#include <signal.h>
#include <sys/prctl.h>
#include <stdexcept>

#ifdef GPERFTOOLS
#include <gperftools/heap-profiler.h>
#include <gperftools/profiler.h>
#endif

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <fstream>
#include <sstream>

#include "modules/bio_base/biograph_dir.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_format/vcf.h"
#include "modules/bio_mapred/flatten_seqset.h"
#include "modules/bio_mapred/merge_flat_seqset.h"
#include "modules/io/autostats.h"
#include "modules/io/command.h"
#include "modules/io/config.h"
#include "modules/io/defaults.h"
#include "modules/io/file_io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/parallel.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/io/version.h"
#include "modules/main/main.h"
#include "modules/mapred/task_mgr.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/align.h"
#include "modules/variants/assemble.h"
#include "modules/variants/calc_coverage.h"
#include "modules/variants/dedup.h"
#include "modules/variants/genotype.h"
#include "modules/variants/normalize.h"
#include "modules/variants/pipeline.h"
#include "modules/variants/ploid_limit.h"
#include "modules/variants/ploidless_vcf_export.h"
#include "modules/variants/ref_map.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/simple_genotype_filter.h"
#include "modules/variants/split_variants.h"
#include "modules/variants/trace_ref.h"

namespace fs = boost::filesystem;
using namespace variants;

class DiscoveryMain : public Main {
 public:
  DiscoveryMain() {
    m_usage =
        "%1% version %2%\n\n"
        "Usage: %1% [OPTIONS] --in <biograph> --ref <reference path> --out "
        "<vcf name> [--sample <sample id>]\n\n"
        "Call variants on a BioGraph.\n";
  }

 protected:
  void add_args() override;
  int run(po::variables_map vars) override;
  const product_version& get_version() override { return biograph_current_version; }
  void warn_memory_cache(const std::string& item) const;
  void write_aligned_csv_assembly_header(writable* out) const;
  void write_aligned_csv_assembly(writable* w, const assemble_options& options, const assembly& a);
  void write_csv_assembly_header(writable* out) const;
  void write_csv_assembly(writable* w, const assemble_options& options, const assembly& a);

 private:
  std::string m_in_biograph;
  std::string m_in_seqset;
  std::string m_readmap_str;
  std::string m_in_readmap;
  std::string m_ref_dir;
  std::string m_assembly_out_file;
  std::string m_aligned_assembly_out_file;
  std::string m_half_aligned_out_file;
  std::string m_vcf_out_file;
  std::string m_bed_file;
  std::string m_chunk_stats_file;

  bool m_force = false;
  bool m_verify_assemble = true;
  bool m_enable_pop_tracer = true;
  bool m_use_bidir_tracer = false;
  bool m_simple_gt = false;
  bool m_rvg_exclude = true;
  bool m_report_long_traces = false;
  unsigned m_min_pop_overlap;
  std::string m_ref_map_file;
  unsigned m_min_overlap;
  float m_min_overlap_pct;
  std::mutex m_csv_mu;
  std::mutex m_aligned_csv_mu;
  int m_max_ploids = 2;

  std::map<std::string, std::string> m_vcf_headers;

  void do_assemble();
  void check_for_terminate();
};

static volatile bool terminate = false;

// handle termination in the main loop.
static void signal_handler(int sig) {
  // One is enough
  signal(sig, SIG_IGN);
  terminate = true;
}

static void update_progress(const float& new_progress) {
  static float prev_progress = 0;
#ifdef GPERFTOOLS
  // Update every percent.
  if (int(new_progress * 100) != int(prev_progress * 100)) {
    ProfilerFlush();
    prev_progress = new_progress;
    print_progress(new_progress);
  }
#endif
  if ((fabs(new_progress - prev_progress) > 0.0001) or (new_progress == 1.0)) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

void DiscoveryMain::check_for_terminate() {
  if (terminate) {
    std::cerr << "\nControl-C detected.\n";
    SPLOG("Control-C detected.");
    m_keep_tmp = true;
    cleanup(false);
    ::exit(1);
  }
}

void DiscoveryMain::add_args() {
  m_general_options.add_options()                                                 //
      ("in", po::value(&m_in_biograph)->required(), "Input BioGraph to process")  //
      ("ref", po::value(&m_ref_dir)->required(), "Reference directory")           //
      ("out", po::value(&m_vcf_out_file)->default_value("-"), "Output VCF file")          //
      ("sample", po::value(&m_readmap_str)->default_value(""),
       "Sample ID (Accession ID, uuid, or coverage file) to process. Only required if the BioGraph "
       "contains multiple samples.")  //
      ("force,f", po::bool_switch(&m_force)->default_value(false),
       "Overwrite existing output file")  //
      ;

  m_variant_options.add_options()  //
      ("bed,regions", po::value(&m_bed_file)->default_value(""),
       "If specified, only call in the regions contained in the given BED file.")  //
      ("min-overlap", po::value(&m_min_overlap_pct)->default_value(.7, ".7"),
       "Minimum overlap required between reads when tracing paths, as a fraction of the read "
       "length (0.5-0.9 recommended)")  //
      ("max-ploids", po::value(&m_max_ploids)->default_value(2),
       "Maximum number of alleles to output")  //
      ;

  m_advanced_options.add_options()  //
      ("assemblies-out", po::value(&m_assembly_out_file)->default_value(""),
       "If specified, assemblies are written to this file in CSV format")  //
      ("aligned-assemblies-out", po::value(&m_aligned_assembly_out_file)->default_value(""),
       "If specified, aligned assemblies are written to this file in CSV format")  //
      ("ref-map", po::value(&m_ref_map_file)->default_value(""),
       "If specified, filename to use to store the reference map between runs.  Warning: No "
       "validity checking is done to make sure this matches the current reference and seqset being "
       "used.")  //
      ("half-aligned-out", po::value(&m_half_aligned_out_file)->default_value(""),
       "If specified, assemblies which are only aligned on one end are written to this file in CSV "
       "format")  //
      ;

  m_secret_options.add_options()  //
      ("verify-assemble", po::value(&m_verify_assemble)->default_value(true),
       "Enable sanity checks when processing assemblies")  //
      ("report-long-traces", po::value(&m_report_long_traces)->default_value(false),
       "Report positions where we spend a long time doing path traces")  //
      ("chunk-stats-out", po::value(&m_chunk_stats_file)->default_value(""),
       "If specified, statistics are written to this file on how long it takes to process each "
       "chunk")  //
      ("enable-pop-tracer", po::value(&m_enable_pop_tracer)->default_value(true),
       "If specified, use a 'pop-front' based tracer in addition to the normal 'push-front-drop' "
       "tracer.")  //
      ("use-bidir-tracer", po::value(&m_use_bidir_tracer)->default_value(false),
       "If specified, use the bidirectional tracer for discovery instead of the older pop and "
       "push tracers.")  //
      ("rvg-exclude", po::value(&m_rvg_exclude)->default_value(true),
       "If specified, exclude low coverage non-structural variants.")  //
      ("verbose-trace-work", po::value(&trace_ref::g_verbose_trace_work)->default_value(false),
       "If true, report in the log whenever traces of regions start or finish")  //
      ("min-pop-overlap",
       po::value(&m_min_pop_overlap)->default_value(assemble_options().min_pop_overlap),
       "Minimum overlap for the pop tracer")  //
      ("simple-gt", po::value(&m_simple_gt)->default_value(false),
       "Attempt simple genotyping and filtering during discovery phase at the expense of "
       "sensitivity")  //
      ;

  m_positional.add("in", 1);
  m_positional.add("ref", 1);
  m_positional.add("out", 1);

  m_options.add(m_general_options).add(m_variant_options);
}

int DiscoveryMain::run(po::variables_map vars) {
  canon_assembly_order::set_default_sort_order(sort_order::OLD_DISCOVER);

  if (fs::exists(m_vcf_out_file) && m_vcf_out_file != "-") {
    if (m_force) {
      fs::remove(m_vcf_out_file);
    } else {
      std::cerr << "Refusing to overwrite '" << m_vcf_out_file << "'. Use -f to override.\n";
      exit(1);
    }
  }
  if (m_assembly_out_file != "" && fs::exists(m_assembly_out_file)) {
    if (m_force) {
      fs::remove(m_assembly_out_file);
    } else {
      std::cerr << "Refusing to overwrite '" << m_assembly_out_file << "'. Use -f to override.\n";
      exit(1);
    }
  }
  if (m_half_aligned_out_file != "" && fs::exists(m_half_aligned_out_file)) {
    if (m_force) {
      fs::remove(m_half_aligned_out_file);
    } else {
      std::cerr << "Refusing to overwrite '" << m_half_aligned_out_file
                << "'. Use -f to override.\n";
      exit(1);
    }
  }

  if (m_stats_file.empty()) {
    m_stats_file = m_in_biograph + "/qc/variants_stats.json";
  }

  if (m_min_overlap_pct < 0.5 || m_min_overlap_pct > 0.9) {
    SPLOG("WARNING: %f overlap is outside of suggested range (0.5, 0.9)", m_min_overlap_pct);
    std::cerr << "WARNING: " << m_min_overlap_pct << " is outside of suggested range (0.5, 0.9)\n";
  }

  initialize_app(m_ref_dir, m_in_biograph + "/qc/variants_log.txt");
  if (m_ref_dir.empty() or not defaults.check_refdir(m_ref_dir)) {
    throw std::runtime_error("Please check your reference directory.");
  }
  // initialize_app() ignores SIGINT, so handle it ourselves.
  signal(SIGINT, signal_handler);

  // Get the seqset and readmap
  biograph_dir bgdir(m_in_biograph, READ_BGDIR);

  m_in_seqset = bgdir.seqset();
  m_in_readmap = bgdir.find_readmap(m_readmap_str);
  m_readmap_str = bgdir.find_readmap_accession(m_readmap_str);

  m_vcf_headers["command-line"] = m_cmdline.c_str();

  // Do the assemble state
  do_assemble();

  // Additional stats
  m_stats.add("command", "variants");
  m_stats.add("version", biograph_current_version.make_string());
  m_stats.add("accession_id", m_readmap_str);
  m_stats.add("reference", m_ref_dir);

  m_stats.save();

  std::cerr << "\n" << m_vcf_out_file << " created.\n";

  return 0;
}

void DiscoveryMain::warn_memory_cache(const std::string& item) const {
  SPLOG("WARNING: %s doesn't seem to be cached in RAM!", item.c_str());
  std::cerr << "WARNING: random access to " << item << " seems slow even after caching it in RAM. "
            << "Does this machine have enough RAM to hold it all?\n";
}

struct variant_stats : public autostats_base {
  DECLARE_AUTOSTATS(variant_stats,               //
                    ((COUNTER, snp))             // total SNPs
                    ((COUNTER, snp_het))         // SNP het
                    ((COUNTER, snp_hom))         // SNP hom
                    ((COUNTER, ts))              // SNP transition (A<->G or C<->T)
                    ((COUNTER, tv))              // SNP transversion (anything else)
                    ((COUNTER, ins))             // total insertions
                    ((COUNTER, ins_het))         // heterozygous ins
                    ((COUNTER, ins_hom))         // homozygous ins
                    ((COUNTER, del))             // total deletions
                    ((COUNTER, del_het))         // heterozygous del
                    ((COUNTER, del_hom))         // homozygous del
                    ((COUNTER, repl))            // total pure replacements
                    ((COUNTER, repl_het))        // pure replacement het
                    ((COUNTER, repl_hom))        // pure replacement hom
                    ((COUNTER, subins))          // total substitution insertions
                    ((COUNTER, subins_het))      // heterozygous subins
                    ((COUNTER, subins_hom))      // homozygous subins
                    ((COUNTER, subdel))          // total substitution deletions
                    ((COUNTER, subdel_het))      // heterozygous subdel
                    ((COUNTER, subdel_hom))      // homozygous subdel
                    ((COUNTER, ins_1))           // ins size 1
                    ((COUNTER, ins_2_9))         // ins 2-9
                    ((COUNTER, ins_10_49))       // ins 10-49
                    ((COUNTER, ins_50_299))      // ins 50-299
                    ((COUNTER, ins_300_999))     // ins 300-999
                    ((COUNTER, ins_1000))        // ins 1000+
                    ((COUNTER, del_1))           // del size 1
                    ((COUNTER, del_2_9))         // del 2-9
                    ((COUNTER, del_10_49))       // del 10-49
                    ((COUNTER, del_50_299))      // del 50-299
                    ((COUNTER, del_300_999))     // del 300-999
                    ((COUNTER, del_1000))        // del 1000+
                    ((COUNTER, repl_2_9))        // repl 2-9
                    ((COUNTER, repl_10_49))      // repl 10-49
                    ((COUNTER, repl_50_299))     // repl 50-299
                    ((COUNTER, repl_300_999))    // repl 300-999
                    ((COUNTER, repl_1000))       // repl 1000+
                    ((COUNTER, subins_1))        // substitution insertion size 1
                    ((COUNTER, subins_2_9))      // subins 2-9
                    ((COUNTER, subins_10_49))    // subins 10-49
                    ((COUNTER, subins_50_299))   // subins 50-299
                    ((COUNTER, subins_300_999))  // subins 300-999
                    ((COUNTER, subins_1000))     // subins 1000+
                    ((COUNTER, subdel_1))        // substitution deletion size 1
                    ((COUNTER, subdel_2_9))      // subdel 2-9
                    ((COUNTER, subdel_10_49))    // subdel 10-49
                    ((COUNTER, subdel_50_299))   // subdel 50-299
                    ((COUNTER, subdel_300_999))  // subdel 300-999
                    ((COUNTER, subdel_1000))     // subdel 1000+
  );
};

class variant_stats_counter : public assemble_pipeline_interface {
 public:
  variant_stats_counter(const assemble_options& options, pipeline_step_t output)
      : m_output(std::move(output)), m_options(options) {}

  void on_assembly(assembly_ptr a) override {
    std::string var_seq = a->seq.as_string();
    std::string ref_seq =
        m_options.scaffold->subscaffold_str(a->left_offset, a->right_offset - a->left_offset);

    size_t var_size = var_seq.size();
    size_t ref_size = ref_seq.size();

    // SNP
    if ((var_size == ref_size) and (var_size == 1)) {
      m_vcf_stats.snp++;
      // genotype
      if (a->strand_count == 1) {
        m_vcf_stats.snp_het++;
      } else {
        m_vcf_stats.snp_hom++;
      }
      // ts/tv. Be sure to count all affected alleles.
      if ((var_seq == "A" and ref_seq == "G") or (var_seq == "G" and ref_seq == "A") or
          (var_seq == "C" and ref_seq == "T") or (var_seq == "T" and ref_seq == "C")) {
        m_vcf_stats.ts += a->strand_count;
      } else {
        m_vcf_stats.tv += a->strand_count;
      }
    }
    // pure replacement
    else if (var_size == ref_size) {
      m_vcf_stats.repl++;
      // genotype
      if (a->strand_count == 1) {
        m_vcf_stats.repl_het++;
      } else {
        m_vcf_stats.repl_hom++;
      }
      // size histogram
      if (var_size >= 1000) {
        m_vcf_stats.repl_1000++;
      } else if (var_size >= 300) {
        m_vcf_stats.repl_300_999++;
      } else if (var_size >= 50) {
        m_vcf_stats.repl_50_299++;
      } else if (var_size >= 10) {
        m_vcf_stats.repl_10_49++;
      } else {
        m_vcf_stats.repl_2_9++;
      }
    }
    // pure insert
    else if (ref_size == 1) {
      m_vcf_stats.ins++;
      // genotype
      if (a->strand_count == 1) {
        m_vcf_stats.ins_het++;
      } else {
        m_vcf_stats.ins_hom++;
      }
      // size histogram (ignore first reference base)
      if (var_size > 1000) {
        m_vcf_stats.ins_1000++;
      } else if (var_size > 300) {
        m_vcf_stats.ins_300_999++;
      } else if (var_size > 50) {
        m_vcf_stats.ins_50_299++;
      } else if (var_size > 10) {
        m_vcf_stats.ins_10_49++;
      } else if (var_size > 2) {
        m_vcf_stats.ins_2_9++;
      } else {
        m_vcf_stats.ins_1++;
      }
    }
    // pure deletion
    else if (var_size == 1) {
      m_vcf_stats.del++;
      // genotype
      if (a->strand_count == 1) {
        m_vcf_stats.del_het++;
      } else {
        m_vcf_stats.del_hom++;
      }
      // size histogram (ignore first reference base)
      if (ref_size > 1000) {
        m_vcf_stats.del_1000++;
      } else if (ref_size > 300) {
        m_vcf_stats.del_300_999++;
      } else if (ref_size > 50) {
        m_vcf_stats.del_50_299++;
      } else if (ref_size > 10) {
        m_vcf_stats.del_10_49++;
      } else if (ref_size > 2) {
        m_vcf_stats.del_2_9++;
      } else {
        m_vcf_stats.del_1++;
      }
    }
    // substitution insertion
    else if (var_size > ref_size) {
      m_vcf_stats.subins++;
      // genotype
      if (a->strand_count == 1) {
        m_vcf_stats.subins_het++;
      } else {
        m_vcf_stats.subins_hom++;
      }
      size_t ins_size = var_size - ref_size;
      // size histogram
      if (ins_size >= 1000) {
        m_vcf_stats.subins_1000++;
      } else if (ins_size >= 300) {
        m_vcf_stats.subins_300_999++;
      } else if (ins_size >= 50) {
        m_vcf_stats.subins_50_299++;
      } else if (ins_size >= 10) {
        m_vcf_stats.subins_10_49++;
      } else if (ins_size >= 2) {
        m_vcf_stats.subins_2_9++;
      } else {
        m_vcf_stats.subins_1++;
      }
    }
    // substitution deletion
    else if (ref_size > var_size) {
      m_vcf_stats.subdel++;
      // genotype
      if (a->strand_count == 1) {
        m_vcf_stats.subdel_het++;
      } else {
        m_vcf_stats.subdel_hom++;
      }
      size_t del_size = ref_size - var_size;
      // size histogram
      if (del_size >= 1000) {
        m_vcf_stats.subdel_1000++;
      } else if (del_size >= 300) {
        m_vcf_stats.subdel_300_999++;
      } else if (del_size >= 50) {
        m_vcf_stats.subdel_50_299++;
      } else if (del_size >= 10) {
        m_vcf_stats.subdel_10_49++;
      } else if (del_size >= 2) {
        m_vcf_stats.subdel_2_9++;
      } else {
        m_vcf_stats.subdel_1++;
      }
    } else {
      throw io_exception("Impossible case in variant_stats_counter");
    }

    m_output->add(std::move(a));
  }
  ~variant_stats_counter() {
    std::lock_guard<std::mutex> l(g_mu);
    global_stats += m_vcf_stats;
  }
  static variant_stats get_global_stats() { return global_stats; }

 private:
  pipeline_step_t m_output;
  variant_stats m_vcf_stats;
  assemble_options m_options;
  static variant_stats global_stats;
  static std::mutex g_mu;
};

variant_stats variant_stats_counter::global_stats;
std::mutex variant_stats_counter::g_mu;

class vcf_pipeline : public scaffold_pipeline_interface {
 public:
  vcf_pipeline(writable* sink) : m_sink(sink) {}

  std::unique_ptr<pipeline_interface> pipeline_for_scaffold(
      const assemble_options& options, const std::string& scaffold_name) override {
    CHECK(options.scaffold);
    CHECK(!options.scaffold_name.empty());
    auto vcf =
        make_unique<ploidless_vcf_export>(options, scaffold_name, [&](const std::string line) {
          static std::mutex mu;
          std::lock_guard<std::mutex> l(mu);
          m_sink->write(line.data(), line.size());
        });
    std::unique_ptr<assemble_pipeline> p = make_unique<assemble_pipeline>(options, std::move(vcf));
    p->add_standard_variants_pipeline();
    p->add_step<variant_stats_counter>(options);
    std::unique_ptr<pipeline_interface> pi = std::move(p);
    return pi;
  }

 private:
  writable* m_sink = nullptr;
};

void DiscoveryMain::write_csv_assembly_header(writable* out) const {
  out->print(
      "scaffold_name,left_offset,right_offset,left_anchor_len,right_anchor_len,"
      "aid,score,min_overlap,ref_seq,seq,generated_by\n");
}

void DiscoveryMain::write_csv_assembly(writable* w, const assemble_options& options,
                                      const assembly& a) {
  std::stringstream out;
  out << options.scaffold_name << "," << a.left_offset << "," << a.right_offset << ","
      << a.left_anchor_len << "," << a.right_anchor_len << "," << a.assembly_id << "," << a.score
      << "," << a.min_overlap << ","
      << options.scaffold->subscaffold_str(a.left_offset, a.right_offset - a.left_offset) << ","
      << a.seq.as_string();
  if (a.tags.empty()) {
    out << ",UNKNOWN";
  } else {
    out << "," << a.tags.to_string_short();
  }
  out << "\n";
  std::string s = out.str();

  std::lock_guard<std::mutex> l(m_csv_mu);
  w->write(s.data(), s.size());
}

void DiscoveryMain::write_aligned_csv_assembly_header(writable* out) const {
  out->print(
      "scaffold_name,left_offset,right_offset,left_anchor_len,right_anchor_len,"
      "aid,score,min_overlap,variants(refrange:varseq:refseq),seq,generated_by\n");
}

void DiscoveryMain::write_aligned_csv_assembly(writable* w, const assemble_options& options,
                                              const assembly& a) {
  std::stringstream out;
  out << options.scaffold_name << "," << a.left_offset << "," << a.right_offset << ","
      << a.left_anchor_len << "," << a.right_anchor_len << "," << a.assembly_id << "," << a.score
      << "," << a.min_overlap << ",";

  bool first_variant = true;
  for (const auto& var : a.aligned_variants) {
    if (first_variant) {
      first_variant = false;
    } else {
      out << ";";
    }
    out << var.left_offset << "-" << var.right_offset << ":" << var.seq.as_string();
    if (options.scaffold) {
      out << ":"
          << options.scaffold->subscaffold_str(var.left_offset, var.right_offset - var.left_offset);
    }
  }
  out << "," << a.seq.as_string() << "," << a.tags.to_string_short() << "\n";
  std::string s = out.str();
  std::lock_guard<std::mutex> l(m_aligned_csv_mu);
  w->write(s.data(), s.size());
}

void DiscoveryMain::do_assemble() {
  // Assembling can sometimes take longer than the default limit of 20
  // minutes.
  Config::set("task_timeout", 3600 * 10);

  m_stats.start_stage("load_seqset");
  SPLOG("Loading seqset: %s", m_in_seqset.c_str());
  spiral_file_options sfopts;
  sfopts.read_into_ram = m_cache_all;
  auto ss = std::make_shared<seqset>(m_in_seqset, sfopts);
  m_vcf_headers["seqset-uuid"] = ss->uuid();
  m_stats.add("uuid", ss->uuid());
  check_for_terminate();

  SPLOG("Caching seqset into RAM");
  std::cerr << "\nLoading biograph\n";

  auto membufs = ss->membufs();
  membufs.cache_in_memory(subprogress(update_progress, 0, 1));
  if (!membufs.is_cached_in_memory()) {
    warn_memory_cache("seqset");
  }
  m_stats.end_stage("load_seqset");
  check_for_terminate();

  SPLOG("Loading readmap: %s", m_in_readmap.c_str());
  readmap rm(ss, m_in_readmap, sfopts);
  if (!rm.has_mate_loop()) {
    throw io_exception("Readmap " + m_in_readmap +
                       " missing mate loop table; upgrade with 'biograph upgrade'");
  }
  check_for_terminate();
  rm.calc_read_len_limits_if_needed();

  SPLOG("Opening reference");
  reference ref("");
  check_for_terminate();

  boost::optional<ref_map> rmap;

  if (m_verify_assemble) {
    assemble_pipeline_interface::global_set_verify_order(true);
  }

  m_stats.start_stage("generate_refmap");
  if (m_ref_map_file.empty()) {
    SPLOG("Generating refmap in memory");
    std::cerr << "\nGenerating refmap\n";
    rmap.emplace(ss.get(), &ref);
    rmap->build(update_progress);
  } else {
    if (!boost::filesystem::exists(m_ref_map_file)) {
      std::string new_refmap = m_ref_map_file + ".new";
      unlink(new_refmap.c_str());
      {
        spiral_file_create_mmap c(new_refmap);
        SPLOG("Generating refmap in %s", new_refmap.c_str());
        std::cerr << "\nGenerating refmap\n";
        ref_map build_rmap(ss.get(), &ref, c.create());
        build_rmap.build(update_progress);
      }
      rename(new_refmap.c_str(), m_ref_map_file.c_str());
    }

    SPLOG("Opening refmap %s", m_ref_map_file.c_str());
    std::cerr << "\nOpening refmap\n";
    spiral_file_open_mmap o(m_ref_map_file, sfopts);
    rmap.emplace(ss.get(), &ref, o.open());
  }
  m_stats.end_stage("generate_refmap");
  check_for_terminate();

  if (!membufs.is_cached_in_memory()) {
    warn_memory_cache("seqset");
  }

  m_min_overlap = static_cast<unsigned int>(m_min_overlap_pct * ss->read_len());
  SPLOG("Using min_overlap of %0.2f * %d = %d", m_min_overlap_pct, ss->read_len(), m_min_overlap);

  assemble_options options;
  options.min_overlap = m_min_overlap;
  options.max_ploids = m_max_ploids;
  options.seqset = ss.get();
  options.readmap = &rm;
  options.ref = &ref;
  options.rmap = &rmap.get();
  options.output_assembly_ids = !m_assembly_out_file.empty();
  options.pop_trace_anchor_drop = m_enable_pop_tracer;
  options.output_ml_features = true;
  options.min_pop_overlap = m_min_pop_overlap;
  options.use_bidir_tracer = m_use_bidir_tracer;
  options.rvg_exclude = m_rvg_exclude;
  options.simple_genotype_filter = m_simple_gt;
  if (m_use_bidir_tracer) {
    options.scaffold_split_size = 400 * 1000;
  }

  std::mutex half_aligned_mu;
  boost::optional<std::ofstream> half_aligned_out;
  if (!m_half_aligned_out_file.empty()) {
    half_aligned_out.emplace(m_half_aligned_out_file);
    (*half_aligned_out) << "scaffold_name,left_anchor,right_anchor,sequence,assembly_id\n";
    options.report_half_aligned_func = [&](const half_aligned_assembly& ha) {
      std::lock_guard<std::mutex> l(half_aligned_mu);
      (*half_aligned_out) << ha.scaffold_name << ",";
      if (ha.right_anchor) {
        (*half_aligned_out) << "," << ha.offset;
      } else {
        (*half_aligned_out) << ha.offset << ",";
      }
      (*half_aligned_out) << "," << ha.seq << "," << ha.assembly_id << "\n";
    };
  }

  if (m_report_long_traces) {
    options.report_long_traces_func = [&](const std::string& scaffold_name, double seconds,
                                          aoffset_t offset, assemble_stats st) {
      static std::mutex mu;
      std::lock_guard<std::mutex> l(mu);
      std::cerr << "Long trace occured on " << scaffold_name << ":" << offset << ", " << seconds
                << " seconds, stats: " << st << "\n";
    };
  }
  std::mutex chunk_stats_mu;
  boost::optional<std::ofstream> chunk_stats_out;
  if (!m_chunk_stats_file.empty()) {
    chunk_stats_out.emplace(m_chunk_stats_file);
    (*chunk_stats_out) << "scaffold_name,start,limit,dir,seconds,stats\n";
    options.report_chunk_stats_func = [&](const std::string& scaffold_name, aoffset_t start,
                                          aoffset_t limit, bool rev_comp, double seconds,
                                          assemble_stats stats) {
      std::lock_guard<std::mutex> l(chunk_stats_mu);
      (*chunk_stats_out) << scaffold_name << "," << start << "," << limit << ","
                         << (rev_comp ? "rev" : "fwd") << "," << seconds << "," << stats << "\n";
      chunk_stats_out->flush();
    };
  }

  std::unique_ptr<file_writer> aligned_assemblies_out;
  if (!m_aligned_assembly_out_file.empty()) {
    aligned_assemblies_out = make_unique<file_writer>(m_aligned_assembly_out_file);
    write_aligned_csv_assembly_header(aligned_assemblies_out.get());
    options.report_aligned_assemblies_func = [&](const assemble_options& options,
                                                 const assembly& a) {
      write_aligned_csv_assembly(aligned_assemblies_out.get(), options, a);
    };
  }

  auto add_ml_features = [&](const assemble_options& options, assembly& a) {
    CHECK(!a.ml_features) << "ml_features should not have already been populated";
    a.ml_features.emplace();
    assembly_ml_features& f = *a.ml_features;
    f.score = a.score;
    f.refspan = a.right_offset - a.left_offset;
    f.lanch = a.left_anchor_len;
    f.ranch = a.right_anchor_len;
    if (f.refspan) {
      int refgc_count = 0;
      for (dna_base b :
           options.scaffold->subscaffold(a.left_offset, a.right_offset - a.left_offset)) {
        if (b == dna_base('G') || b == dna_base('C')) {
          ++refgc_count;
        }
      }
      f.refgc = float(refgc_count) / f.refspan;
    } else {
      f.refgc = 0;
    }
    int asmlen = a.seq.size();
    if (asmlen) {
      int altgc_count = 0;
      for (dna_base b : a.seq) {
        if (b == dna_base('G') || b == dna_base('C')) {
          ++altgc_count;
        }
      }
      f.altgc = float(altgc_count) / asmlen;
    } else {
      f.altgc = 0;
    }
    f.alt_seq = a.seq;
  };

  std::unique_ptr<file_writer> assemblies_out;
  if (!m_assembly_out_file.empty()) {
    assemblies_out = make_unique<file_writer>(m_assembly_out_file);
    write_csv_assembly_header(assemblies_out.get());
    options.report_discovered_assemblies_func = [&](const assemble_options& options, assembly& a) {
      write_csv_assembly(assemblies_out.get(), options, a);
      add_ml_features(options, a);
    };
  } else {
    options.report_discovered_assemblies_func = add_ml_features;
  }

  std::unique_ptr<file_writer> vcf_out;
  if (m_vcf_out_file == "-")
  {
  	vcf_out = make_unique<file_writer>("/dev/stdout");
  } else {
  	vcf_out = make_unique<file_writer>(m_vcf_out_file);
  }
  auto header = ploidless_vcf_export::header(options, m_vcf_headers, m_readmap_str);
  vcf_out->write(header.data(), header.size());
  auto p = make_unique<vcf_pipeline>(vcf_out.get());

  trace_ref t(options, p.get());
#ifdef GPERFTOOLS
  ProfilerStart("/scratch/biograph_variants.prof");
//  HeapProfilerStart("/scratch/biograph_variants_heap.prof");
#endif

  if (m_bed_file.empty()) {
    SPLOG("Assembling whole reference");
    t.add_entire_reference();
  } else {
    std::cerr << "\nAssembling extents in BED file: " << m_bed_file << "\n";
    SPLOG("Assembling extents in BED file %s", m_bed_file.c_str());

    file_reader bed(m_bed_file);
    std::string line;

    std::vector<std::string> bed_lines;
    while (bed.readline(line, 1000)) {
      std::vector<std::string> fields;
      boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);
      CHECK_GE(fields.size(), 3) << "bad BED line: '" << line << "'";

      std::string scaffold = fields[0];
      size_t start = atol(fields[1].c_str());
      size_t limit = atol(fields[2].c_str());

      t.add_scaffold_range(scaffold, start, limit);
    }
  }
  std::cerr << "\nAssembling...\n";
  m_stats.start_stage("assemble");
  auto st = t.assemble(update_progress);
  m_stats.end_stage("assemble");

  variant_stats vstats = variant_stats_counter::get_global_stats();

  js::Object report;
  for (const auto& stat : vstats.value_map()) {
    report.push_back(js::Pair(stat.first, stat.second));
  }
  report.push_back(js::Pair(
      "total", vstats.snp + vstats.del + vstats.ins + vstats.repl + vstats.subins + vstats.subdel));
  report.push_back(js::Pair("total_het", vstats.snp_het + vstats.del_het + vstats.ins_het +
                                             vstats.repl_het + vstats.subins_het +
                                             vstats.subdel_het));
  report.push_back(js::Pair("total_hom", vstats.snp_hom + vstats.del_hom + vstats.ins_hom +
                                             vstats.repl_hom + vstats.subins_hom +
                                             vstats.subdel_hom));

  m_stats.add("calls", report);

#ifdef GPERFTOOLS
  ProfilerStop();
// HeapProfilerStop();
#endif

  std::stringstream msg;
  msg << st;
  SPLOG("%s", msg.str().c_str());
}

std::unique_ptr<Main> discovery_main() { return std::unique_ptr<Main>(new DiscoveryMain); }


class AssembleMain : public Main {
 public:
  AssembleMain() {}

  int run(po::variables_map vars) {
    std::cerr << "The 'variants' command has been retired. Please use 'discovery' instead.\n";
    return 1;
  }
};

std::unique_ptr<Main> assemble_main() { return std::unique_ptr<Main>(new AssembleMain); }

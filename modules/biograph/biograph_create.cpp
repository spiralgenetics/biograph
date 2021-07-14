#include <signal.h>
#include <sys/prctl.h>
#include <fstream>
#include <stdexcept>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/regex.hpp>

#include "base/base.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/dna_base_set.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_format/corrected_reads.h"
#include "modules/bio_format/vcf.h"
#include "modules/bio_mapred/kmerize_bf.h"
#include "modules/bio_mapred/make_readmap.h"
#include "modules/bio_mapred/mem_seqset.h"
#include "modules/bio_mapred/read_correction.h"
#include "modules/build_seqset/builder.h"
#include "modules/build_seqset/correct_reads.h"
#include "modules/build_seqset/expand.h"
#include "modules/build_seqset/kmer_counter.h"
#include "modules/build_seqset/part_repo.h"
#include "modules/build_seqset/read_importer.h"
#include "modules/io/config.h"
#include "modules/io/defaults.h"
#include "modules/io/digest.h"
#include "modules/io/file_io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/runtime_stats.h"
#include "modules/io/stopwatch.h"
#include "modules/io/track_mem.h"
#include "modules/io/uuid.h"
#include "modules/io/utils.h"
#include "modules/io/zip.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/manifest_parallel.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/resource_manager.h"
#include "modules/mapred/task_mgr.h"
#include "modules/pipeline/paired_merger.h"

#include "modules/bio_base/biograph_dir.h"
#include "modules/main/main.h"

#include "tools/malloc_select.h"

#if TCMALLOC
# include "gperftools/malloc_hook.h"
# include "gperftools/malloc_extension.h"
#else
# if !defined(__clang__)
#  warning "BioGraph should be compiled with TCMALLOC"
# endif
#endif

using namespace boost::algorithm;
namespace fs = boost::filesystem;

static volatile bool terminate = false;

class SEQSETMain : public Main {
 public:
  SEQSETMain();

 private:
  int run(po::variables_map vars) override;
  void add_args() override;
  po::variables_map get_args();
  void rm_files(const std::string& pattern);
  void do_read_correction(std::unique_ptr<kmer_set> ks, const manifest& uncorrected_reads,
                          manifest& corrected_reads, const read_correction_params& rcp);
  void make_seqset(const manifest& corrected_reads);
  void do_readmap(manifest input_manifest);

  const product_version& get_version() override { return biograph_current_version; }

  class read_importer_state : public parallel_local {
   public:
    struct params {
      std::string tmp_dir;
      build_seqset::kmer_counter* kmer_counter = nullptr;
      bool allow_long_reads;
      std::string tmp_encoding;
      double sample_reads = 0;
      std::mutex* output_mu = nullptr;
      manifest* output_manifest = nullptr;
    };

    using init_type = params;

    read_importer_state() = delete;
    read_importer_state(const read_importer_state&) = delete;
    read_importer_state(params p) : m_counter(*p.kmer_counter), m_params(p), m_uuid(make_uuid()) {}

    void open() {
      CHECK(!m_sink);
      CHECK(!m_sink_f);
      CHECK(!m_raw_sink_f);
      static std::atomic<size_t> g_file_idx{0};
      size_t file_idx = g_file_idx.fetch_add(1);
      CHECK(!m_params.tmp_dir.empty());
      m_path = m_params.tmp_dir + "/all_reads_" + m_uuid + "_" + std::to_string(file_idx);
      m_raw_sink_f = make_unique<file_writer>(m_path);
      m_sink_f = make_encoder(m_params.tmp_encoding, *m_raw_sink_f);
      m_sink = make_unique<kv_writer>(*m_sink_f);
      m_size_written = 0;
      m_records_written = 0;
    }

    void process(const std::vector<std::pair<read_id, unaligned_reads>>& reads) {
      for (const auto& read : reads) {
        if (!m_sink) {
          open();
        }

        if (m_params.sample_reads) {
          m_sample_accum += m_params.sample_reads;
          if (m_sample_accum > 1) {
            m_sample_accum -= 1;
          } else {
            continue;
          }
        }
        for (const auto& r : read.second) {
          if (!m_params.allow_long_reads &&
              r.sequence.size() > std::numeric_limits<uint8_t>::max()) {
            throw(io_exception(printstring(
                "Encountered read of length %ld, which is larger than the maximum read length %d",
                r.sequence.size(), std::numeric_limits<uint8_t>::max())));
          }
          m_counter.add(r.sequence);
        }
        std::string key = msgpack_serialize(read.first);
        std::string value = msgpack_serialize(read.second);
        m_sink->write(key, value);
        m_size_written += key.size() + value.size();
        ++m_records_written;
        if (m_size_written > k_target_size) {
          flush_chunk();
        }
      }
    }

    void flush_chunk() {
      if (!m_sink) {
        return;
      }

      CHECK(m_sink);
      m_sink->close();
      m_sink.reset();
      CHECK(m_sink_f);
      m_sink_f->close();
      m_sink_f.reset();
      CHECK(m_raw_sink_f);
      m_raw_sink_f->close();
      m_raw_sink_f.reset();

      file_info fi(m_path, m_size_written, m_records_written);
      m_local_manifest.add(fi, 0);
      m_path.clear();
    }

    void flush() override {
      flush_chunk();
      std::lock_guard<std::mutex> l(*m_params.output_mu);
      m_params.output_manifest->add(m_local_manifest, true /* unsorted */);
    }

   private:
    std::string m_path;
    std::string m_replaced_sequence;
    std::unique_ptr<kv_writer> m_sink;
    std::unique_ptr<writable> m_raw_sink_f;
    std::unique_ptr<writable> m_sink_f;
    manifest m_local_manifest;
    static constexpr size_t k_target_size = 128 * 1024 * 1024;  // 128 MB
    size_t m_size_written = 0;
    size_t m_records_written = 0;
    build_seqset::kmer_counter::prob_pass_processor m_counter;
    params m_params;
    double m_sample_accum = 0.5;
    std::string m_uuid;
  };

  std::vector<std::string> m_in_reads;
  std::vector<std::string> m_in_pairs;
  std::string m_in_format;
  std::string m_tmp_encoding;
  std::string m_kmer_size;
  std::string m_min_kmer_count;
  std::string m_out;
  std::string m_ref_dir;
  std::string m_min_corrected_reads;
  std::string m_warn_corrected_reads;
  std::string m_accession_id;
  std::string m_readmap_sha;
  std::string m_trim_after_portion;
  std::string m_max_corrections;
  std::string m_min_good_run;
  std::string m_overrep_thresh;
  std::string m_sys_err_thresh;
  std::string m_rnd_err_thresh;
  std::string m_sample_reads;
  std::string m_cut_reads;
  std::string m_dump_kmers;

  bool m_force;
  bool m_allow_long_reads = false;
  bool m_fastq_interleaved = false;
  bool m_got_paired = false;

  size_t m_read_count = 0;
  int m_partition_depth = 0;

  biograph_dir m_bgdir;

  // Progress updater that also checks for termination.
  progress_handler_t m_update_progress;

  std::unique_ptr<build_seqset::part_counts> m_part_counts;
};

constexpr size_t SEQSETMain::read_importer_state::k_target_size;

SEQSETMain::SEQSETMain() {
  m_update_progress = [&](double new_progress) -> void {
    static float prev_progress = 0;
    if (fabs(new_progress - prev_progress) > 0.0001) {
      prev_progress = new_progress;
      print_progress(new_progress);
    }
    if (terminate) {
      std::cerr << "\nControl-C detected.\n";
      SPLOG("Control-C detected.");
      m_keep_tmp = true;
      cleanup(false);
      ::exit(1);
    }
  };

  m_usage =
      "%1% version %2%\n\n"
      "Usage: %1% [OPTIONS] --reads <file> --ref <refdir> --out <biograph> "
      "[--pair <fastq pairs>] [...]\n\n"
      "Convert reads to BioGraph format.";
}

void SEQSETMain::add_args() {
  m_general_options.add_options()("out", po::value(&m_out)->required(),
                                  "Output BioGraph name (.bg)")          //
      ("ref", po::value(&m_ref_dir)->required(), "Reference directory")  //
      ("reads,in", po::value(&m_in_reads)->required(),
       "Input file to process (fastq, bam, cram. Use - for STDIN)")  //
      ("format", po::value(&m_in_format)->default_value("auto"),
       "Input file format when using STDIN (fastq, bam, cram)")  //
      ("interleaved", po::bool_switch(&m_fastq_interleaved)->default_value(false),
       "Input reads are interleaved (fastq only)")                                              //
      ("pair", po::value(&m_in_pairs), "Second input file containing read pairs (fastq only)")  //
      ("id", po::value(&m_accession_id)->default_value(""),
       "Optional accession ID for this sample")  //
      ("force,f", po::bool_switch(&m_force)->default_value(false),
       "Overwrite existing BioGraph")  //
      ;

  m_kmer_options.add_options()  //
      ("min-kmer-count", po::value(&m_min_kmer_count)->default_value("5"),
       "The integer minimum kmer count. Reads with kmers less abundant "
       "than this will be corrected or dropped. (min 1)")  //
      ("kmer-size", po::value(&m_kmer_size)->default_value("30"),
       "The size of kmers to use for kmer generation.")  //
      ("trim-after-portion", po::value(&m_trim_after_portion)->default_value("0.7"),
       "Trim the end of reads until they pass read correction, down to a minimum of the given "
       "portion of the read length. 1 = no automatic trimming.")  //
      ;

  m_correction_options.add_options()  //
      ("max-corrections", po::value(&m_max_corrections)->default_value("8"),
       "Correct up to the specified number of bases.")  //
      ("min-good-run", po::value(&m_min_good_run)->default_value("2"),
       "Minimum number of good bases between corrections.")  //
      ("min-reads", po::value(&m_min_corrected_reads)->default_value("0.4"),
       "Minimum fraction of reads that must survive read correction")  //
      ("warn-reads", po::value(&m_warn_corrected_reads)->default_value("0.7"),
       "Warn when this fraction of reads does not survive read correction")  //
      ;

  m_advanced_options.add_options()  //
      ("allow-long-reads", po::bool_switch(&m_allow_long_reads)->default_value(false),
       "Allow reads longer than 255 bases (EXPERIMENTAL)")  //
      ("tmp-encoding", po::value(&m_tmp_encoding)->default_value("gzip1"),
       "Encoding to use for temporary files; using \"gzip\" here will use more CPU time and less "
       "I/O.  \"null\" means store temporary files uncompressed.  \"gzip1\" specifies a "
       "compression level of 1, which is faster than default but doesn't compress as well.")  //
      ;
  // Options such as --max-mem:
  m_general_options.add(track_mem_program_options());

  m_secret_options.add_options()  //
      ("overrep-threshold", po::value(&m_overrep_thresh)->default_value("0"),
       "If non-zero, the number of times a kmer must occur before attempting overrepresentation "
       "filtering.  If zero, overrepresentation filtering is disabled.")  //
      ("sys-err-thresh", po::value(&m_sys_err_thresh)->default_value("0.1"),
       "Systematic error threshold for overrep filtering")  //
      ("rnd-err-thresh", po::value(&m_rnd_err_thresh)->default_value("0.005"),
       "Rnd error threshold for overrep filtering")  //
      ("sample-reads", po::value(&m_sample_reads)->default_value("0"),
       "If non-zero, sample this portion of the input reads; other reads are ignored.  For "
       "instance, to sample 1 in 4 reads, specify --sample-reads=0.25")  //
      ("cut-reads", po::value(&m_cut_reads)->default_value(""),
       "If present, drop all bases in each read that aren't in this range.  For example, "
       "--cut-reads=10-100 will only use the 10th through the 100th base of each read, "
       "resulting in read lengths of up to 90 bases.  (The resultant read lengths could "
       "be shorter if the input reads are less than 100 bases long to start with.)")  //
      ("dump-kmers", po::value(&m_dump_kmers)->default_value(""),
       "If present, output an unsorted kmer list to the given file for use in the kmer set "
       "benchmark")  //
      ;

  m_options.add(m_general_options);
  m_advanced_options.add(m_kmer_options).add(m_correction_options);

  m_positional.add("in", 1);
  m_positional.add("ref", 1);
  m_positional.add("out", 1);
}

static size_t validate_param(const std::string& param, const std::string& value,
                             const std::vector<size_t> range = {}) {
  size_t num_val;
  try {
    num_val = std::stoull(value);
  } catch (const std::exception& ex) {
    throw std::runtime_error(boost::str(boost::format("%s must specify an integer") % param));
  }
  if ((range.size() >= 1) && (num_val < range[0])) {
    throw std::runtime_error(
        boost::str(boost::format("%s must specify an integer >= %d") % param % range[0]));
  }
  if ((range.size() == 2) && (num_val > range[1])) {
    throw std::runtime_error(
        boost::str(boost::format("%s must specify an integer <= %d") % param % range[1]));
  }
  return num_val;
}

static float validate_float_param(const std::string& param, const std::string& value,
                                  const std::vector<float> range = {}) {
  float num_val;
  try {
    num_val = std::stof(value);
  } catch (const std::exception& ex) {
    throw std::runtime_error(
        boost::str(boost::format("%s must specify a floating point number") % param));
  }
  if ((range.size() >= 1) && (num_val < range[0])) {
    throw std::runtime_error(boost::str(
        boost::format("%s must specify a floating point number >= %f") % param % range[0]));
  }
  if ((range.size() == 2) && (num_val > range[1])) {
    throw std::runtime_error(boost::str(
        boost::format("%s must specify a floating point number <= %f") % param % range[1]));
  }
  return num_val;
}

static std::pair<unsigned, unsigned> validate_cut_param(const std::string& param,
                                                        const std::string& value) {
  if (value.empty()) {
    return std::make_pair(0, 0);
  }
  auto it = value.find('-');
  if (it == std::string::npos) {
    throw(std::runtime_error(
        boost::str(boost::format("%s must specify a range separated by a dash") % param)));
  }
  std::string start_str = value.substr(0, it);
  unsigned start_val;
  try {
    start_val = std::stoul(start_str);
  } catch (const std::exception& ex) {
    throw(std::runtime_error(boost::str(
        boost::format("%s must specify a numerical range; couldn't parse %s as a number") % param %
        start_str)));
  }
  std::string end_str = value.substr(it + 1);
  unsigned end_val;
  try {
    end_val = std::stoul(end_str);
  } catch (const std::exception& ex) {
    throw(std::runtime_error(boost::str(
        boost::format("%s must specify a numerical range; couldn't parse %s as a number") % param %
        end_str)));
  }

  if (end_val <= start_val) {
    throw(std::runtime_error(
        boost::str(boost::format("%s must specify a nonzero range; %d must be less than %d") %
                   param % start_val % end_val)));
  }

  return std::make_pair(start_val, end_val);
}

// termination is handled in the main loop.
static void signal_handler(int sig) {
  // One is enough
  signal(sig, SIG_IGN);
  terminate = true;
}

// Delete files under m_tmp_dir that match pattern (simple substring match)
void SEQSETMain::rm_files(const std::string& pattern) {
  fs::directory_iterator end_itr;
  for (fs::directory_iterator i(m_tmp_dir); i != end_itr; i++) {
    if (i->path().filename().string().find(pattern) != std::string::npos) {
      fs::remove(i->path());
    }
  }
}

int SEQSETMain::run(po::variables_map vars) {
#if TCMALLOC
# if NDEBUG
  constexpr bool k_ndebug = true;
# else
  constexpr bool k_ndebug = false;
# endif

  if (!k_ndebug) {
    // Having this compiled whether or not we're in NDEBUG helps
    // ensure that we don't accidentally leave out TCMALLOC, which
    // could cause severe performance problems running biograph.

    MallocHook::AddNewHook(track_mem::get_malloc_new_hook());
  }

  track_mem::set_reset_stats_hook([]() {
    auto inst = MallocExtension::instance();

    // Return all memory from this stage to OS immediately.
    size_t unmapped_before = 0;
    CHECK(inst->GetNumericProperty("tcmalloc.pageheap_unmapped_bytes", &unmapped_before))
        << "Unable to get tmcalloc statistics";

    auto ms = stopwatch([inst]() { inst->ReleaseFreeMemory(); }).count();

    size_t unmapped_after = 0;
    CHECK(inst->GetNumericProperty("tcmalloc.pageheap_unmapped_bytes", &unmapped_after))
        << "Unable to get tmcalloc statistics";

    constexpr ssize_t k_mb = 1024 * 1024;
    ssize_t unmapped_diff = ssize_t(unmapped_after) - ssize_t(unmapped_before);
    if (unmapped_diff > k_mb || ms > 1000) {
      SPLOG("Returned %ld MB free memory to OS in %ld ms (total unmapped = %ld MB)",
            unmapped_diff / k_mb, ms, unmapped_after / k_mb);
    }
  });
#endif

  // Reference is required for many stages
  if (not defaults.check_refdir(m_ref_dir)) {
    throw std::runtime_error("Please check your reference directory.");
  }

  size_t minimum_mem = 32;
  size_t configured_max_gb = get_maximum_mem_bytes() / 1024 / 1024 / 1024;

  if (configured_max_gb < minimum_mem) {
    std::cerr << "WARNING: Configured memory limit of " << configured_max_gb
              << " GiB is less than recommended " << minimum_mem << " GiB\n";
  }

  size_t kmer_size = validate_param("--kmer-size", m_kmer_size, {16, 32});
  size_t min_kmer_count = validate_param("--min-kmer-count", m_min_kmer_count, {1, 10000000});
  float min_corrected_reads =
      validate_float_param("min-reads", m_min_corrected_reads, {0.0, 1.0});
  float warn_corrected_reads =
      validate_float_param("warn-reads", m_warn_corrected_reads, {0.0, 1.0});
  float trim_after_portion =
      validate_float_param("trim-after-portion", m_trim_after_portion, {0.0, 1.0});
  unsigned max_corrections = validate_param("max-corrections", m_max_corrections, {0, 32});
  unsigned min_good_run = validate_param("min-good-run", m_min_good_run, {0, 64});

  unsigned overrep_thresh = validate_param("overrep-threshold", m_overrep_thresh, {0, 10000000});
  float sys_err_thresh = validate_float_param("sys-err-threshold", m_sys_err_thresh, {0.0, 1.0});
  float rnd_err_thresh = validate_float_param("rnd-err-threshold", m_rnd_err_thresh, {0.0, 1.0});
  float sample_reads = validate_float_param("sample-reads", m_sample_reads, {0.0, 1.0});
  std::pair<unsigned, unsigned> cut_reads = validate_cut_param("cut-reads", m_cut_reads);

  std::set<std::string> formats = {"bam", "cram", "fastq", "auto"};
  if (formats.find(m_in_format) == formats.end()) {
    std::cerr << "Invalid input format '" << m_in_format << "'\n";
    ::exit(1);
  }

  if (!m_in_pairs.empty() && m_in_pairs.size() != m_in_reads.size()) {
    std::cerr << "If pair files are present, there must be the same number of them as read files.\n";
    ::exit(1);
  }
  // We don't have forcing on and there's an existing file
  if (!m_force && biograph_dir::force_check(m_out)) {
    std::cerr << "Refusing to overwrite '" + m_out + "'. Use --force to override.\n";
    ::exit(1);
  }

  m_bgdir = biograph_dir(m_out, CREATE_BGDIR);

  if (m_accession_id.empty()) {
    m_accession_id = fs::canonical(fs::path(m_out)).stem().string();
  }

  if (m_stats_file.empty()) {
    m_stats_file = m_out + "/qc/create_stats.json";
  }

  // Initialize and kick off the daemons
  initialize_app(m_ref_dir, m_out + "/qc/create_log.txt");

  // Now set up the custom handler
  signal(SIGINT, signal_handler);
  signal(SIGTERM, signal_handler);

  manifest reads;
  std::vector<manifest> kmers_and_hist;
  manifest corrected;
  manifest coverage;

  reference ref("", m_ref_dir);

  m_stats.start_stage("import");
  std::cerr << "Importing reads\n";

  build_seqset::count_kmer_options kmer_opts;
  kmer_opts.kmer_size = kmer_size;
  kmer_opts.max_memory_bytes = get_maximum_mem_bytes();
  kmer_opts.min_count = min_kmer_count;
  kmer_opts.progress = subprogress(m_update_progress, 0.0, 0.05);

  // Don't waste memory and storage making a probabilistic table
  // that's going to be less than 1% full.
  kmer_opts.max_prob_table_entries = ref.size() * 100;

  build_seqset::kmer_counter counter(kmer_opts);
  read_importer_state::params import_params;
  import_params.tmp_dir = m_tmp_dir;
  import_params.kmer_counter = &counter;
  import_params.allow_long_reads = m_allow_long_reads;
  import_params.tmp_encoding = m_tmp_encoding;
  import_params.sample_reads = sample_reads;
  import_params.output_manifest = &reads;
  std::mutex output_manifest_mu;
  import_params.output_mu = &output_manifest_mu;
  build_seqset::read_importer<read_importer_state> importer(
      import_params, subprogress(m_update_progress, 0.05, 1.0));
  if (cut_reads.second) {
    importer.set_cut_region(cut_reads.first, cut_reads.second);
  }

  if (!m_in_pairs.empty()) {
    CHECK_EQ(m_in_reads.size(), m_in_pairs.size());
  }

  for (size_t i = 0; i != m_in_reads.size(); ++i) {
    std::string in_reads = m_in_reads[i];
    std::string in_pairs;
    if (!m_in_pairs.empty()) {
      in_pairs = m_in_pairs[i];
    }

    std::string in_format = m_in_format;

    if (in_reads == "-") {
      in_reads = "/dev/stdin";
    }

    if (in_format == "auto") {
      if (m_fastq_interleaved) {
        std::cerr << "--interleaved specified. Assuming fastq format.\n";
        in_format = "fastq";
      } else if (!in_pairs.empty()) {
        std::cerr << "--pair specified. Assuming fastq format.\n";
        in_format = "fastq";
      } else if (ends_with(in_reads, ".bam")) {
        in_format = "bam";
      } else if (ends_with(in_reads, ".cram")) {
        in_format = "cram";
      } else if (ends_with(in_reads, ".fq") or ends_with(in_reads, ".fq.gz") or
                 ends_with(in_reads, ".fastq") or ends_with(in_reads, ".fastq.gz")) {
        in_format = "fastq";
      } else if (in_reads == "/dev/stdin") {
        in_format = "bam";
        std::cerr << "Streaming reads from STDIN with no --format specified. Assuming bam or cram.\n";
      } else {
        std::cerr << "Cannot determine the input file type of " << in_reads << ".\n";
        std::cerr << "Input file does not end in .bam .cram .fq .fastq .fq.gz or .fastq.gz.\n";
        std::cerr << "Please specify --format.\n";
        ::exit(1);
      }
    }

    if (in_format == "fastq") {
      importer.queue_fastq(in_reads, in_pairs, m_fastq_interleaved);
    } else {
      if (!in_pairs.empty()) {
        std::cerr << "Non-fastq input " << in_reads
                  << " may not containing separate pairing data in " << in_pairs << "\n";
        ::exit(1);
      }
      importer.queue_bam(in_reads, m_ref_dir);
    }
  }

  SPLOG("Initializing kmer counter");
  counter.start_prob_pass();
  SPLOG("Importing reads");

  m_read_count = importer.import();
  m_got_paired = importer.got_paired();
  reads.metadata().set(meta::ns::readonly, "paired", m_got_paired);
  reads.set_encoding(m_tmp_encoding);

  if (!m_in_pairs.empty() && !m_got_paired) {
    throw std::runtime_error("Pair files specified but no pairs were successfully imported");
  }

  if (m_read_count == 0) {
    bool has_cram = false;
    for (const auto& in_reads : m_in_reads) {
      if (ends_with(in_reads, ".cram")) {
        has_cram = true;
      }
    }
    if (has_cram) {
      std::cerr << "\nCheck that the reference " << m_ref_dir << " matches the input cram.\n";
    }
    throw std::runtime_error("\nNo reads were imported, exiting.");
  }
  print_progress(1.0);

  std::cerr << "\nTotal reads imported: " << m_read_count << std::endl;
  SPLOG("%lu reads imported", m_read_count);

  size_t m_imported_count = 0;
  for (const auto& i : reads) {
    m_imported_count += i.num_records;
  }
  if (sample_reads) {
    std::cerr << "After sampling and pair association, we have " << m_read_count << " reads\n";
    SPLOG("%lu reads present after sampling and pair association", m_imported_count);
  } else {
    SPLOG("%lu reads present after pair association", m_imported_count);
  }
  m_stats.end_stage("import");

  m_stats.start_stage("kmerization");
  std::cerr << "\nRunning kmerization\n";

  // Limit kmerization to 2 threads small datasets due to contention
  if (m_read_count < 50 * 1000 * 1000 and get_thread_count() > 2) {
    SPLOG("Small dataset. Limiting kmerization to 2 threads.");
    set_thread_count("2");
  }

  counter.set_progress_handler(subprogress(m_update_progress, 0.0, 0.85));
  counter.close_prob_pass();
  set_thread_count(m_requested_threads);
  SPLOG("done close_prob_pass");

  kmerize_bf_params kbf;
  kbf.kmer_size = kmer_size;
  kbf.error_rate = 0.05;
  kbf.reference = "";
  kbf.memory_bound = get_maximum_mem_bytes() / 1024 / 1024 / 1024;
  kbf.num_threads = m_num_threads;
  kbf.min_count = min_kmer_count;
  kbf.overrep = overrep_thresh;
  kbf.sys_err_thresh = sys_err_thresh;
  kbf.rnd_err_thresh = rnd_err_thresh;
  kbf.dump_kmers_file = m_dump_kmers;

  std::unique_ptr<kmer_set> ks;
  std::tie(ks, kmers_and_hist) =
      run_kmerize_subtask(kbf, reads, &counter, subprogress(m_update_progress, 0.85, 1.0));
  CHECK(ks);

  // Kmerization
  CHECK_EQ(3, kmers_and_hist.size());
  kmers_and_hist[0].update_metadata(reads);

  fs::copy_file(m_tmp_dir + "/kmer_quality_report.html", m_out + "/qc/kmer_quality_report.html",
                fs::copy_option::overwrite_if_exists);
  for (std::string fr : get_kmer_filter_result_types()) {
    if (fs::exists(m_tmp_dir + "/kmer_quality_report-" + fr + ".html")) {
      fs::copy_file(m_tmp_dir + "/kmer_quality_report-" + fr + ".html",
                    m_out + "/qc/kmer_quality_report-" + fr + ".html",
                    fs::copy_option::overwrite_if_exists);
    }
  }

  print_progress(1.0);
  m_stats.end_stage("kmerization");

  // Read correction
  m_stats.start_stage("read_correction");
  std::cerr << "\nCorrecting reads\n";

  if (m_read_count < 10 * 1000 * 1000) {
    m_partition_depth = 2;
  } else if (m_read_count < 100 * 1000 * 1000) {
    m_partition_depth = 3;
  } else {  // m_read_count >= 100 * 1000 * 1000
    m_partition_depth = 4;
  }
  SPLOG("Using a partition depth of %d (%d partitions)", m_partition_depth,
        1 << (2 * m_partition_depth));

  read_correction_params rcp;
  rcp.min_kmer_score = min_kmer_count;
  rcp.skip_snps = false;
  rcp.exact = (max_corrections == 0);
  rcp.trim_after_portion = trim_after_portion;
  rcp.frc_max_corrections = max_corrections;
  rcp.frc_min_good_run = min_good_run;

  do_read_correction(std::move(ks), reads, corrected, rcp);
  corrected.update_metadata(kmers_and_hist[0]);

  size_t num_corrected_bases =
      corrected.metadata().get<uint64_t>(meta::ns::readonly, "corrected_read_bases", 0);
  float cov_estimate = num_corrected_bases * 1. / ref.size();
  SPLOG("%0.2fx estimated corrected coverage", cov_estimate);

  size_t num_corrected_reads =
      corrected.metadata().get(meta::ns::readonly, "corrected_read_count", 0);

  float corrected_pct = num_corrected_reads * 1. / m_read_count;
  if (corrected_pct < min_corrected_reads) {
    std::string msg =
        boost::str(boost::format("Fewer than %2.0f%% of reads (set by "
                                 "--min-reads) were kept "
                                 "after correction (%lu / %lu remain). "
                                 "Cannot continue.") %
                   (min_corrected_reads * 100.0) % num_corrected_reads % m_read_count);
    SPLOG("%s", msg.c_str());
    throw std::runtime_error(msg);
  }
  if (corrected_pct < warn_corrected_reads) {
    std::string msg =
        boost::str(boost::format("Warning: Fewer than %2.0f%% of reads (set by "
                                 "--warn-reads) survived correction "
                                 "(%lu / %lu remain)") %
                   (warn_corrected_reads * 100.0) % num_corrected_reads % m_read_count);
    SPLOG("%s", msg.c_str());
    std::cerr << msg << "\n";
  } else {
    SPLOG("%lu / %lu reads survived read correction.", num_corrected_reads, m_read_count);
  }

  if (not m_keep_tmp) {
    SPLOG("Deleting kmers");
    rm_files("kmerize_");
  }

  print_progress(1.0);
  m_stats.end_stage("read_correction");

  m_stats.start_stage("make_seqset");
  make_seqset(corrected);
  m_stats.end_stage("make_seqset");

  m_stats.start_stage("make_readmap");
  do_readmap(corrected);
  m_stats.end_stage("make_readmap");

  m_stats.start_stage("metadata");

  // Save metadata
  biograph_metadata meta = m_bgdir.get_metadata();

  meta.accession_id = m_accession_id;

  samples_t samples = {{m_accession_id, m_readmap_sha}};
  meta.samples = samples;

  m_bgdir.set_metadata(meta);
  m_bgdir.save_metadata();

  m_stats.add("command", "create");
  m_stats.add("version", biograph_current_version.make_string());
  m_stats.add("accession_id", m_accession_id);
  m_stats.add("reference", m_ref_dir);
  m_stats.add("imported_reads", m_read_count);
  m_stats.add("coverage", cov_estimate);
  m_stats.add("corrected_reads", num_corrected_reads);
  m_stats.add("corrected_bases", num_corrected_bases);
  m_stats.add("avg_bases_per_read", num_corrected_bases * 1. / num_corrected_reads);
  m_stats.add("corrected_pct", corrected_pct);
  m_stats.add("uuid", m_bgdir.biograph_id());

  m_stats.save();
  m_stats.end_stage("metadata");

  std::cerr << "\n" << m_out << " created.\n";

  return 0;
}

void SEQSETMain::do_readmap(manifest input_manifest) {
  std::cerr << "\nCalculating coverage...\n";

  spiral_file_options sfopts;
  seqset_file the_seqset_file(m_bgdir.seqset(), sfopts.with_read_into_ram(m_cache_all));
  std::string readmap_path = m_bgdir.readmap("tmp");
  make_readmap::do_make(readmap_path, the_seqset_file, input_manifest, m_got_paired,
                        the_seqset_file.get_seqset().max_read_len(), m_update_progress);

  m_readmap_sha = sha1sum(fs::path(readmap_path));
  fs::rename(readmap_path, m_bgdir.readmap(m_readmap_sha));

  print_progress(1.0);
}

std::unique_ptr<Main> seqset_main() { return std::unique_ptr<Main>(new SEQSETMain); }

void SEQSETMain::do_read_correction(std::unique_ptr<kmer_set> ks, const manifest& uncorrected,
                                    manifest& corrected, const read_correction_params& rcp) {
  CHECK(ks);
  SPLOG("Fast creation enabled");
  boost::optional<build_seqset::part_repo> entries;
  entries.emplace(m_partition_depth, m_tmp_dir + "/seq_ref-", m_tmp_dir + "/seq_repo");

  reference ref("", m_ref_dir);
  size_t ref_size = ref.size();
  SPLOG("Found %ld bases of reference", ref_size);
  entries->add_initial_repo(dna_slice(ref.get_dna(0), ref.size()));

  build_seqset::correct_reads cr(*entries, *ks, rcp);
  cr.add_initial_repo(subprogress(m_update_progress, 0, 0.1));
  SPLOG("Correcting reads...");
  entries->open_write_pass("initial");
  std::vector<file_info> file_infos;
  for (auto the_file_info : uncorrected) {
    file_infos.push_back(the_file_info);
  }
  size_t corrected_read_count = 0;
  size_t corrected_read_bases = 0;
  std::mutex m_cr_out_mu;
  parallel_for(  //
      0, file_infos.size(),
      [&](size_t idx) {
        const auto& fi = file_infos[idx];
        auto file_reader_ptr = fi.file.read();
        auto file_decoder = make_decoder(uncorrected.get_encoding(), *file_reader_ptr);
        kv_reader file_kv_reader(*file_decoder);
        read_id key;
        unaligned_reads value;

        manifest local_corrected;
        output_stream_params cr_osp;
        cr_osp.encoding = m_tmp_encoding;
        std::unique_ptr<kv_sink> cr_sink(
            cr_osp.build(CONF_S(path_bulkdata), "corrected_reads", local_corrected));
        size_t local_corrected_count = 0;
        size_t local_corrected_bases = 0;

        corrected_reads c;
        while (file_kv_reader.read_msgpack(key, value)) {
          for (const auto& r : value) {
            c.emplace_back();
            if (cr.correct(r, c.back())) {
              local_corrected_bases += c.back().corrected.size();
            } else {
              c.pop_back();
            }
          }
          if (!c.empty()) {
            cr_sink->write_msgpack(key.pair_name, c);
            local_corrected_count += c.size();
            c.clear();
          }
        }
        cr_sink->close();
        cr_sink.reset();

        // Remove uncorrected reads; we don't need them anymore.
        path(fi.file).remove();

        std::lock_guard<std::mutex> l(m_cr_out_mu);
        corrected.add(local_corrected, true);
        corrected_read_count += local_corrected_count;
        corrected_read_bases += local_corrected_bases;
      },
      subprogress(m_update_progress, 0.1, 1));

  SPLOG("Generated %ld corrected reads, %ld bases (avg %.2f bases/read)", corrected_read_count,
        corrected_read_bases, corrected_read_bases * 1. / corrected_read_count);
  corrected.metadata().set(meta::ns::readonly, "corrected_read_count", corrected_read_count);
  corrected.metadata().set(meta::ns::readonly, "corrected_read_bases", corrected_read_bases);

  m_part_counts = entries->release_part_counts("initial");
  CHECK(m_part_counts);
}

void SEQSETMain::make_seqset(const manifest& corrected) {
  SPLOG("Fast creation enabled");
  boost::optional<build_seqset::part_repo> entries;
  entries.emplace(m_partition_depth, m_tmp_dir + "/seq_ref-", m_tmp_dir + "/seq_repo");
  std::cerr << "\nGenerating BioGraph\n";
  entries->flush();
  entries->reset_part_counts("initial", std::move(m_part_counts));
  {
    equal_subprogress expand_progress(subprogress(m_update_progress, 0, 0.8), 4);
    SPLOG("Expanding");
    build_seqset::expander expand(*entries, m_keep_tmp);
    expand.sort_and_dedup("", "initial", "init_sorted", "", 0, 0, expand_progress[0]);
    expand.expand("init_sorted", "init_expanded", 7, 255, expand_progress[1]);
    expand.sort_and_dedup("init_sorted", "init_expanded", "pass2_sorted", "pass2_expanded", 1, 6,
                          expand_progress[2]);
    expand.sort_and_dedup("pass2_sorted", "pass2_expanded", "complete", "", 0, 0,
                          expand_progress[3]);
  }
  SPLOG("Building seqset");
  track_mem::reset_stats();
  build_seqset::builder b;
  b.build_chunks(*entries, "complete", m_keep_tmp, subprogress(m_update_progress, 0.8, 0.9));

  track_mem::reset_stats();
  if (not m_keep_tmp) {
    entries->partitions("complete", false /* don't need pushed iterators */,
                        true /* delete on close */);
    rm_files("seq_repo");
  }
  entries.reset();
  {
    spiral_file_create_mmap c(m_out + "/seqset");
    b.make_seqset(c.create(), subprogress(m_update_progress, 0.9, 1));
  }

  print_progress(1.0);
}

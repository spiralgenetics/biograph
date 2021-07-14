#include <signal.h>
#include <sys/prctl.h>
#include <stdexcept>

#include <limits.h>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <vector>

#include "modules/bio_base/make_mergemap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_base/seqset_mergemap.h"
#include "modules/bio_base/seqset_merger.h"
#include "modules/bio_mapred/flatten_seqset.h"
#include "modules/bio_mapred/make_readmap.h"
#include "modules/bio_mapred/merge_flat_seqset.h"
#include "modules/io/config.h"
#include "modules/io/digest.h"
#include "modules/io/file_io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/io/version.h"

#include "modules/bio_base/biograph_dir.h"
#include "modules/main/main.h"

namespace js = json_spirit;
namespace fs = boost::filesystem;

class MergeSEQSETMain : public Main {
 public:
  MergeSEQSETMain() {
    m_usage =
        "%1% version %2%\n\n"
        "Usage: %1% [OPTIONS] --out <merged biograph> --in <source biograph> <source biograph> "
        "[...]\n\n"
        "Merge BioGraphs. Produces a single merged BioGraph with coverage data for every\n"
        "sample in the input BioGraphs.\n";
  }

 protected:
  void add_args() override;
  int run(po::variables_map vars) override;
  const product_version& get_version() override { return biograph_current_version; }

 private:
  // tmp_path("/foo/bar.txt", ".xyz") -> "/tmp/spiral_XyZZy/bar.xyz"
  std::string tmp_path(const std::string& in_file, const std::string& extension);
  void do_merge();
  void warn_memory_cache(const std::string& item) const;

  po::options_description m_seqpath_options{"Assembly Options"};

  std::vector<std::string> m_in_files;
  std::string m_out;
  std::string m_accession_id;
  bool m_force;
  bool m_use_full_ids;
  // size_t m_min_overlap;

  biograph_dir m_out_bgdir;
  std::vector<biograph_dir> m_in_bgdirs;
  samples_t m_samples;
};

std::string temp_directory;

// termination is handled in the main loop.
static void signal_handler(int sig) {
  // One is enough
  signal(sig, SIG_IGN);
  std::cerr << "\nControl-C detected.\n";
  std::cerr << "Temp directory retained in " << temp_directory << std::endl;
  ::exit(1);
}

static void update_progress(const float& new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

void MergeSEQSETMain::add_args() {
  m_general_options.add_options()                                                           //
      ("out", po::value(&m_out)->required(), "Output merged BioGraph")                      //
      ("in", po::value(&m_in_files)->required()->multitoken(), "Input biographs to merge")  //
      ("id", po::value(&m_accession_id)->default_value(""),
       "Optional accession ID for the merged BioGraph")  //
      ("force,f", po::bool_switch(&m_force)->default_value(false),
       "Overwrite existing BioGraph")  //
      ;

  m_positional.add("out", 1);
  m_positional.add("in", -1);

  m_options.add(m_general_options);
}

int MergeSEQSETMain::run(po::variables_map vars) {
  std::set<std::string> in_bg_ids;
  std::set<std::string> in_bg_accession_ids;
  std::set<std::string> in_sample_ids;

  // pre-flight check: iterate over all input biographs
  for (auto& in_file : m_in_files) {
    biograph_dir bgdir(in_file, READ_BGDIR);
    if (not bgdir.is_valid()) {
      throw std::runtime_error("Cannot open '" + in_file + "': invalid BioGraph.");
    }

    // Check that all biographs are unique
    if (in_bg_ids.find(bgdir.biograph_id()) != in_bg_ids.end()) {
      std::cerr << "Duplicate BioGraph ID for '" << in_file << "', skipping.\n";
      continue;
    }

    if (in_bg_accession_ids.find(bgdir.accession_id()) != in_bg_accession_ids.end()) {
      std::cerr << "Duplicate Accession ID '" << bgdir.accession_id() << "' for '" << in_file
                << "', skipping.\n";
      continue;
    }

    // Check that all sample ids are unique
    if (bgdir.samples().empty()) {
      // auto-assign IDs HERE
      throw std::runtime_error("No sample metadata found for '" + in_file + "'. Cannot continue.");
    }

    in_bg_ids.insert(bgdir.biograph_id());
    in_bg_accession_ids.insert(bgdir.accession_id());
    m_in_bgdirs.push_back(bgdir);

    // TODO: glob readmaps, check if 1:1
    // readmap but no metadata? accession_id = "readmap's biograph accession_id: sample_1 ..."

    // If any sample accession is reused, use fully qualified accession IDs
    if (not m_use_full_ids) {
      for (auto sample : bgdir.samples()) {
        if (in_sample_ids.find(sample.first) != in_sample_ids.end()) {
          m_use_full_ids = true;
          break;
        }
        in_sample_ids.insert(sample.first);
      }
    }
  }

  if (m_in_bgdirs.size() < 2) {
    throw std::runtime_error("Merge requires two or more unique BioGraphs.");
  }

  if (!m_force && biograph_dir::force_check(m_out)) {
    std::cerr << "Refusing to overwrite '" + m_out + "'. Use --force to override.\n";
    exit(1);
  }
  m_out_bgdir = biograph_dir(m_out, CREATE_BGDIR);

  if (m_stats_file.empty()) {
    m_stats_file = m_out + "/qc/merge_stats.json";
  }

  initialize_app("", m_out + "/qc/merge_log.txt");
  temp_directory = m_tmp_dir;

  // initialize_app() ignores SIGINT, so handle it ourselves.
  signal(SIGINT, signal_handler);

  do_merge();

  return 0;
}

std::string MergeSEQSETMain::tmp_path(const std::string& in_file, const std::string& extension) {
  return fs::path(fs::path(m_tmp_dir) / fs::path(in_file).stem()).string() + extension;
}

void MergeSEQSETMain::warn_memory_cache(const std::string& item) const {
  SPLOG("WARNING: %s doesn't seem to be cached in RAM!", item.c_str());
  std::cerr << "WARNING: random access to " << item
            << " seems slow even after caching it in RAM.  Does this machine "
               "have enough RAM to hold it all?\n";
}

void MergeSEQSETMain::do_merge() {
  m_stats.start_stage("make_flats");
  for (auto& in_bgdir : m_in_bgdirs) {
    std::string flat_path = tmp_path(in_bgdir.biograph_id(), ".flat");

    std::cerr << in_bgdir.path() << std::endl;
    SPLOG("Building flat seqset for %s", in_bgdir.path().c_str());
    std::unique_ptr<seqset_file> ss_f(new seqset_file(in_bgdir.seqset()));

    SPLOG("Caching %s into RAM", in_bgdir.path().c_str());
    auto membufs = ss_f->membufs();
    membufs.cache_in_memory(subprogress(update_progress, 0, 0.05));
    if (!membufs.is_cached_in_memory()) {
      warn_memory_cache(in_bgdir.path());
    }

    SPLOG("Creating spiral file");
    std::unique_ptr<spiral_file_create_mmap> sp_mmap(new spiral_file_create_mmap(flat_path));

    SPLOG("Creating flat output");
    std::unique_ptr<seqset_flat_builder> flat(new seqset_flat_builder(&ss_f->get_seqset()));

    SPLOG("Building flat");
    flat->build(sp_mmap->create(), subprogress(update_progress, 0.05, 1.0));

    SPLOG("Flat build complete");
    flat.reset();
    sp_mmap->close();

    std::cerr << std::endl;
  }

  SPLOG("Opening flats");
  std::vector<std::unique_ptr<seqset_file>> seqsets;
  std::vector<const seqset_flat*> flat_ptrs;

  for (auto& in_bgdir : m_in_bgdirs) {
    std::string flat_path = tmp_path(in_bgdir.biograph_id(), ".flat");
    seqsets.emplace_back(new seqset_file(in_bgdir.seqset()));
    spiral_file_open_mmap sp_mmap(flat_path);
    flat_ptrs.emplace_back(new seqset_flat(sp_mmap.open(), &seqsets.back()->get_seqset()));
  }
  m_stats.end_stage("make_flats");

  m_stats.start_stage("make_mergemaps");
  spiral_file_create_mmap create_merge(m_out_bgdir.seqset());

  SPLOG("Building mergemaps");
  make_mergemap mm_make(flat_ptrs);
  std::cerr << "Creating merge maps" << std::endl;
  mm_make.build(subprogress(update_progress, 0.0, 0.95));

  SPLOG("%lu entries in resultant merge; writing mergemaps", mm_make.total_merged_entries());

  unsigned input_index = 0;
  float prev = 0.95;
  float inc = 0.05 / m_in_files.size();
  for (auto& in_bgdir : m_in_bgdirs) {
    std::string in_path = tmp_path(in_bgdir.biograph_id(), ".mergemap");
    SPLOG("Building mergemap %s", in_path.c_str());
    spiral_file_create_mmap sp_mmap(in_path);
    seqset_mergemap_builder build_mergemap(sp_mmap.create(),
                                           seqsets[input_index]->get_seqset().uuid(),
                                           create_merge.uuid(), mm_make.total_merged_entries());
    mm_make.fill_mergemap(input_index++, &build_mergemap,
                          subprogress(update_progress, prev, prev + inc));
    prev = prev + inc;

    sp_mmap.close();
  }
  print_progress(1.0);
  std::cerr << std::endl;

  SPLOG("Opening mergemaps");
  std::vector<const seqset_mergemap*> mergemap_ptrs;
  for (auto& in_bgdir : m_in_bgdirs) {
    spiral_file_open_mmap sp_mmap(tmp_path(in_bgdir.biograph_id(), ".mergemap"));
    SPLOG("adding %s", tmp_path(in_bgdir.biograph_id(), ".mergemap").c_str());
    mergemap_ptrs.emplace_back(new seqset_mergemap(sp_mmap.open()));
  }
  m_stats.end_stage("make_mergemaps");

  m_stats.start_stage("final_merge");
  seqset_merger merger(flat_ptrs, mergemap_ptrs);

  std::cerr << "Generating merged BioGraph" << std::endl;
  merger.build(create_merge.create(), update_progress);
  print_progress(1.0);
  std::cerr << std::endl;
  create_merge.close();

  // Put your toys away. They're quite large.
  for (auto& flat_ptr : flat_ptrs) {
    delete (flat_ptr);
  }
  for (auto& mergemap_ptr : mergemap_ptrs) {
    delete (mergemap_ptr);
  }
  m_stats.end_stage("final_merge");

  m_stats.start_stage("create_readmaps");
  for (auto& in_bgdir : m_in_bgdirs) {
    for (const auto& sample : in_bgdir.samples()) {
      SPLOG("Migrating %s:%s", in_bgdir.biograph_id().c_str(), sample.second.c_str());
      std::cerr << "Coverage: " << in_bgdir.accession_id() << " (" << sample.first << ")"
                << std::endl;

      std::string readmap_path = in_bgdir.readmap(sample.second);
      std::string mergemap_path = tmp_path(in_bgdir.biograph_id(), ".mergemap");
      std::string output_path = m_out_bgdir.readmap("tmp");
      fs::remove(output_path);

      SPLOG("Opening original readmap");
      auto old_readmap = readmap::open_anonymous_readmap(readmap_path);
      spiral_file_open_mmap o(mergemap_path);

      SPLOG("Opening mergemap");
      seqset_mergemap mergemap(o.open());

      SPLOG("Opening new readmap");
      spiral_file_create_mmap new_readmap(output_path);

      SPLOG("Everything opened, starting migration");
      make_readmap::fast_migrate(*old_readmap, mergemap, new_readmap.create(), update_progress);

      std::string sha = sha1sum(fs::path(output_path));
      SPLOG("Rename tmp readmap to %s", m_out_bgdir.readmap(sha).c_str());
      fs::rename(output_path, m_out_bgdir.readmap(sha));

      if (m_use_full_ids) {
        m_samples[in_bgdir.accession_id() + ":" + sample.first] = sha;
      } else {
        m_samples[sample.first] = sha;
      }

      print_progress(1.0);
      std::cerr << std::endl;
    }
  }
  m_stats.end_stage("create_readmaps");

  m_stats.start_stage("metadata");
  std::vector<std::string> command_history;
  biograph_metadata meta = m_out_bgdir.get_metadata();
  if (m_accession_id.empty()) {
    for (auto& in_bgdir : m_in_bgdirs) {
      if (m_accession_id.empty()) {
        m_accession_id = in_bgdir.accession_id();
      } else {
        m_accession_id = m_accession_id + "+" + in_bgdir.accession_id();
      }
    }
  }

  // Command History and logs/QC
  for (auto& in_bgdir : m_in_bgdirs) {
    spiral_file_open_mmap sf(in_bgdir.seqset());

    command_history.push_back(sf.file_info().command_line_str());
    for (auto const& cmd : in_bgdir.get_metadata().command_history) {
      command_history.push_back(cmd);
    }

    fs::path qc_in(in_bgdir.path() / fs::path("qc/create_log.txt"));
    if (fs::exists(qc_in)) {
      fs::copy_file(qc_in, fs::path(m_out_bgdir.path() /
                                    fs::path("qc/" + in_bgdir.accession_id() + "_create_log.txt")));
    }
    fs::path kq_in(in_bgdir.path() / fs::path("qc/kmer_quality_report.html"));
    if (fs::exists(kq_in)) {
      fs::copy_file(kq_in, fs::path(m_out_bgdir.path() / fs::path("qc/" + in_bgdir.accession_id() +
                                                                  "_kmer_quality_report.html")));
    }

    // I'm worried this is going to lose tracking information or it'll have overwriting or something
    fs::directory_iterator end_itr;
    for (fs::directory_iterator i(fs::path(in_bgdir.path() / fs::path("qc"))); i != end_itr; i++) {
      if (i->path().filename().string().find("_log.txt") != std::string::npos) {
        fs::path dest(m_out_bgdir.path() / fs::path("qc/") /
                      fs::path(in_bgdir.accession_id() + "_" + i->path().filename().string()));
        if (!fs::exists(dest)) {
          fs::copy_file(i->path(), dest);
        }
      }
      if (i->path().filename().string().find("_kmer_quality_report.html") != std::string::npos) {
        fs::path dest(m_out_bgdir.path() / fs::path("qc/") /
                      fs::path(in_bgdir.accession_id() + "_" + i->path().filename().string()));
        if (!fs::exists(dest)) {
          fs::copy_file(i->path(), dest);
        }
      }
    }
  }

  meta.accession_id = m_accession_id;
  meta.samples = m_samples;
  meta.command_history = command_history;

  m_out_bgdir.set_metadata(meta);
  m_out_bgdir.save_metadata();

  m_stats.add("command", "merge");
  m_stats.add("version", biograph_current_version.make_string());
  m_stats.add("accession_id", m_accession_id);
  m_stats.add("samples", m_samples.size());
  m_stats.add("uuid", m_out_bgdir.biograph_id());

  m_stats.save();

  m_stats.end_stage("metadata");

  std::cerr << std::endl << m_out << " created." << std::endl;
}

std::unique_ptr<Main> merge_seqset_main() { return std::unique_ptr<Main>(new MergeSEQSETMain); }

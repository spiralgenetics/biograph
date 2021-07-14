#include <signal.h>
#include <sys/prctl.h>
#include <stdexcept>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include "modules/bio_base/biograph_dir.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_mapred/make_readmap.h"
#include "modules/io/config.h"
#include "modules/io/digest.h"
#include "modules/io/file_io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/version.h"

#include "modules/main/main.h"

static void update_progress(const float& new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

class UpgradeReadmapMain : public Main {
 public:
  UpgradeReadmapMain() {
    m_usage =
        "%1% version %2%\n\n"
        "Usage: %1% [OPTIONS] --in <target biograph>\n\n"
        "Upgrades all readmaps in a biograph from v2 (mate pairs) to v3 (mate "
        "loops)\n";
    ;
  }

 protected:
  void add_args() override;
  int run(po::variables_map vars) override;
  const product_version& get_version() override { return biograph_current_version; }

 private:
  std::string m_bgdir;
};

// anchored handles termination in the main loop.
static void signal_handler(int sig) {
  // One is enough
  signal(sig, SIG_IGN);
  std::cout << "\nControl-C detected.\n";
  ::exit(1);
}

void UpgradeReadmapMain::add_args() {
  m_general_options.add_options()("in", po::value(&m_bgdir)->required(), "Target biograph");
  m_options.add(m_general_options);

  m_positional.add("in", 1);
}

int UpgradeReadmapMain::run(po::variables_map vars) {
  initialize_app("", m_bgdir + "/qc/upgrade_readmap_log.txt");

  // initialize_app() ignores SIGINT, so handle it ourselves.
  signal(SIGINT, signal_handler);

  std::cout << "Opening biograph" << std::endl;

  SPLOG("Opening biograph %s", m_bgdir.c_str());
  biograph_dir bgdir(m_bgdir, READ_BGDIR);
  if (!bgdir.is_valid()) {
    throw io_exception(m_bgdir + " is not a valid biograph");
  }
  SPLOG("Opening seqset %s", bgdir.seqset().c_str());
  std::shared_ptr<seqset> ss_f = std::make_shared<seqset>(bgdir.seqset());

  SPLOG("Caching %s into RAM", bgdir.path().c_str());
  ss_f->membufs().cache_in_memory(subprogress(update_progress, 0, 1.0));
  print_progress(1.0);

  boost::optional<seqset_flat> flat;

  samples_t out_samples;

  for (const auto& sample : bgdir.samples()) {
    SPLOG("Migrating %s:%s", bgdir.biograph_id().c_str(), sample.second.c_str());
    std::string readmap_path = bgdir.readmap(sample.second);
    readmap old_readmap(ss_f, readmap_path);

    if (old_readmap.has_mate_loop()) {
      std::cout << "\n" << sample.first << " is already upgraded" << std::endl;
      SPLOG("%s already has mate loop enabled", sample.second.c_str());
      out_samples[sample.first] = sample.second;
      continue;
    }

    if (!old_readmap.has_mate_loop()) {
      if (!flat) {
        std::cout << "\nFlattening biograph" << std::endl;
        std::string flat_path =
            fs::path(fs::path(m_tmp_dir) / fs::path(bgdir.biograph_id()).stem()).string() + ".flat";
        {
          std::unique_ptr<spiral_file_create_mmap> sp_mmap(new spiral_file_create_mmap(flat_path));
          SPLOG("Creating flat output");
          std::unique_ptr<seqset_flat_builder> flat(new seqset_flat_builder(&ss_f->get_seqset()));

          SPLOG("Building flat");
          flat->build(sp_mmap->create(), update_progress);

          SPLOG("Flat build complete");
          flat.reset();
          sp_mmap->close();
        }
        print_progress(1.0);
        SPLOG("Opening flat path %s", flat_path.c_str());
        spiral_file_open_mmap sp_mmap(flat_path);
        flat.emplace(sp_mmap.open(), &ss_f->get_seqset());
      }
    }

    if (!old_readmap.has_pairing_data()) {
      std::cout << "\n" << sample.first << " has no pairing data, skipping" << std::endl;
      SPLOG("%s has no pairing data", sample.second.c_str());
      out_samples[sample.first] = sample.second;
      continue;
    }

    std::cout << "\nMigrating " << sample.first << std::endl;
    SPLOG("%ld readmap entries to migrate", old_readmap.size());
    std::string output_path = readmap_path + ".upgraded";
    make_readmap::upgrade(old_readmap, *ss_f, output_path,
                          [&flat](uint64_t seqset_id, unsigned len) -> dna_sequence {
                            CHECK(flat)
                                << "Lookup shouldn't be called if we're not upgrading mate loops";
                            dna_slice seq = flat->get(seqset_id);
                            return dna_sequence(seq.begin(), seq.begin() + len);
                          },
                          subprogress(update_progress, 0.7, 1));

    std::string sha = sha1sum(fs::path(output_path));
    SPLOG("Rename tmp readmap to %s", bgdir.readmap(sha).c_str());
    fs::rename(output_path, bgdir.readmap(sha));

    out_samples[sample.first] = sha;

    fs::remove(fs::path(readmap_path));
    print_progress(1.0);
  }

  SPLOG("updating metadata");

  biograph_metadata meta = bgdir.get_metadata();

  meta.samples = out_samples;
  meta.version = biograph_current_version.make_string();
  meta.command_history.push_back(m_cmdline);

  bgdir.set_metadata(meta);
  bgdir.save_metadata();

  std::cout << "\nUpgrade complete." << std::endl;

  return 0;
}

std::unique_ptr<Main> upgrade_readmap_main() {
  return std::unique_ptr<Main>(new UpgradeReadmapMain);
}

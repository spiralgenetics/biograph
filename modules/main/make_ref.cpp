#include "modules/main/main.h"

#include <signal.h>
#include <sys/prctl.h>
#include <sys/stat.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <stdexcept>

#include "modules/bio_base/flat_ref.h"
#include "modules/io/command.h"
#include "modules/io/config.h"
#include "modules/io/defaults.h"
#include "modules/io/file_io.h"
#include "modules/io/utils.h"
#include "modules/io/zip.h"
#include "modules/mapred/task_mgr.h"
#include "modules/pipeline/build_reference_task.h"

class MakeRefMain : public Main {
 public:
  MakeRefMain() {
    m_usage =
        "%1% version %2%\n\n"
        "Usage: %1% [OPTIONS] --in <fasta> --refdir <ref dir>\n\n"
        "Prepare a fasta for use with BioGraph. The specified reference directory will be\n"
        "created and will contain the new reference and database files.\n";
  }

 protected:
  void add_args() override;
  int run(po::variables_map vars) override;
  const product_version& get_version() override { return biograph_current_version; }

 private:
  std::string m_in_file;
  std::string m_ref_dir;
  std::string m_out_file;
  size_t m_min_n_run = 50;
  bool m_force;
  void check_for_terminate();
  void update_progress(const float& new_progress);
};

namespace fs = boost::filesystem;

static volatile bool terminate = false;

// call this to see if it's time to terminate
void MakeRefMain::check_for_terminate() {
  if (terminate) {
    std::cout << "\nControl-C detected.\n";
    SPLOG("Control-C detected.");
    m_keep_tmp = true;
    cleanup(false);
    ::exit(1);
  }
}

void MakeRefMain::add_args() {
  m_general_options.add_options()("in", po::value(&m_in_file)->required(),
                                  "Input reference fasta or fasta.gz")(
      "refdir", po::value(&m_ref_dir)->required(), "Output reference directory")(
      "force,f", po::bool_switch(&m_force)->default_value(false), "Overwrite existing reference")(
      "min-n-run", po::value(&m_min_n_run)->default_value(50),
      "Any runs of 'N's smaller than this long are replaced with the preceeding base");

  m_positional.add("in", 1);
  m_positional.add("refdir", 1);

  m_options.add(m_general_options);
}

// handle termination in the main loop.
void makeref_signal_handler(int sig) {
  // One is enough
  signal(sig, SIG_IGN);
  terminate = true;
}

// helper to only update progress when the delta is > 0.01%
void MakeRefMain::update_progress(const float& new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

int MakeRefMain::run(po::variables_map vars) {
  if (fs::exists(m_ref_dir)) {
    if (not m_force) {
      throw std::runtime_error("Refusing to overwrite '" + m_ref_dir +
                               "'. Use --force to override.\n");
    }
    if (not fs::exists(fs::path(m_ref_dir) / fs::path("source.fasta"))) {
      throw std::runtime_error(
          m_ref_dir +
          " is not a BioGraph reference. Remove it manually or specify a different location.");
    }
    fs::remove_all(fs::path(m_ref_dir));
  }
  fs::create_directories(m_ref_dir);

  std::string dest = fs::path(fs::path(m_ref_dir) / fs::path(defaults.original_fasta)).string();
  if (boost::algorithm::ends_with(m_in_file, ".gz")) {
    std::cout << "Unzipping source fasta" << std::endl;
    file_reader raw_in_first(m_in_file);
    zip_reader raw_in(raw_in_first);
    file_writer out(dest);

    /*write blocksize for efficiency*/
    struct stat fi;
    stat(m_ref_dir.c_str(), &fi);
    size_t blksz = fi.st_blksize;
    char blk[blksz];

    while (int len = raw_in.read(blk, blksz)) {
      out.write(blk, len);
    }

    out.close();
    raw_in.close();
    raw_in_first.close();
  } else {
    std::cout << "Preparing source fasta" << std::endl;
    fs::copy_file(m_in_file, dest);
  }
  fs::permissions(dest, fs::others_read | fs::owner_read | fs::group_read);

  // Initialize and kick off the daemons
  initialize_app(m_ref_dir);
  launch_daemons();

  // Now set up the custom handler
  signal(SIGINT, makeref_signal_handler);
  signal(SIGTERM, makeref_signal_handler);

  std::cout << "Building reference" << std::endl;
  std::unique_ptr<build_reference_task> tat = make_unique<build_reference_task>(m_ref_dir, "");
  tat->m_min_n_run = m_min_n_run;

  task_mgr tm(new_taskdb_couch());
  std::string id = tm.add_job(CONF_S(path_bulkdata), std::move(tat), "make_ref");
  int job_state = 0;
  float cur_progress = 0.0;

  print_progress(0.0);
  while (job_state == 0) {
    // Time to quit?
    check_for_terminate();
    job_state = tm.state(id);
    cur_progress = tm.get_progress(id);
    update_progress(cur_progress);
    sleep(1);
  }
  print_progress(1.0);
  std::cout << "\n";
  if (job_state != 1) {
    throw io_exception("Reference build could not be completed.");
  }

  std::cout << "Results saved to " << m_ref_dir << "\n";

  return 0;
}

std::unique_ptr<Main> make_ref_main() { return std::unique_ptr<Main>(new MakeRefMain); }

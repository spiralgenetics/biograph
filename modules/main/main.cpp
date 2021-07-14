#include <signal.h>
#include <sys/prctl.h>
#include <sys/utsname.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cstring>
#include <functional>
#include <map>
#include <thread>

#include <boost/algorithm/string/join.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "Poco/Net/ServerSocket.h"
#include "Poco/Net/TCPServer.h"

#include "base/command_line.h"
#include "modules/io/config.h"
#include "modules/io/log.h"
#include "modules/io/parallel.h"
#include "modules/main/main.h"
#include "modules/mapred/path.h"
#include "modules/web/httpserver.h"

std::unique_ptr<Main> manager_main();
std::unique_ptr<Main> taskdb_main();

using Poco::Net::ServerSocket;
using Poco::Net::TCPServer;

namespace fs = boost::filesystem;
namespace pt = boost::posix_time;

namespace {

std::string getproctitle() {
  char buf[2048];
  prctl(PR_GET_NAME, buf);
  return std::string(buf);
}

#if !defined(GPERFTOOLS)
static void log_dumper(int pipe_fd, int log_fd) {
  size_t buffer_size = 65536;
  char buffer[buffer_size];
  while (true) {
    int r = read(pipe_fd, buffer, buffer_size);
    if (r < 0) {
      perror("Unable to read from log pipe");
      return;
    }
    if (r == 0) break;
    r = write(log_fd, buffer, r);
    if (r < 0) {
      perror("Unable to write to log file");
      return;
    }
  }
}

int setup_logger(const std::string& log_file, const bool write_debug) {
  // Make sure the target log dir exists before writing there
  std::string log_path = fs::path(log_file).parent_path().string();
  if (not fs::exists(log_path)) {
    fs::create_directories(log_path);
  }

  int log_fd = open(log_file.c_str(), O_WRONLY | O_APPEND | O_CREAT, 0777);
  if (log_fd < 0) {
    perror(std::string("Unable to open log file " + log_file).c_str());
    exit(1);
  }
  int pipefds[2];
  if (pipe(pipefds) < 0) {
    perror("Unable to create log pipe");
    exit(1);
  }
  int pid_logger = fork();
  if (pid_logger < 0) {
    perror("Unable to fork logger process");
    exit(1);
  }
  if (pid_logger == 0) {
    setproctitle("biograph_logger");
    close(pipefds[1]);
    log_dumper(pipefds[0], log_fd);
    exit(1);
  }
  close(pipefds[0]);
  log_init(getproctitle().c_str(), pipefds[1], write_debug);
  log_build_stamp();

  return pid_logger;
}
#endif

}  // namespace

po::variables_map Main::parse_args(int argc, char* argv[]) {
  std::vector<std::string> arglist(argv, argv + argc);
  m_cmdline = boost::algorithm::join(arglist, " ");

  m_general_options.add_options()                     //
      ("help,h", "Display this help message")         //
      ("help-all", "Show help for advanced options")  //
      ("tmp", po::value(&m_tmp_dir)->default_value(""),
       "Basepath to temporary space. Defaults to a random directory under /tmp/")  //
      ;
  m_advanced_options.add_options()  //
      ("keep-tmp", po::bool_switch(&m_keep_tmp)->default_value(false),
       "Retain temp directory after completion (for debugging)")  //
      ("threads", po::value(&m_requested_threads)->default_value("auto"),
       "Number of concurrent worker threads")                                                    //
      ("debug", po::bool_switch(&m_debug_log)->default_value(false), "Turn on verbose logging")  //
      ("cache", po::bool_switch(&m_cache_all)->default_value(false),
       "Attempt to cache as much as possible in RAM")                                         //
      ("stats", po::value(&m_stats_file)->default_value(""), "Save JSON stats to this file")  //
#if GPERFTOOLS
      ("cpuprofile-dir", po::value(&m_cpuprofile_dir)->default_value(""),
       "Save CPU profiles for each stage to this directory")  //
#endif
      ;

  po::variables_map vars;
  try {
    add_args();

    m_all_options.add(m_options).add(m_advanced_options).add(m_secret_options);

    po::store(
        po::command_line_parser(argc, argv).positional(m_positional).options(m_all_options).run(),
        vars);

    if (vars.count("help-all")) {
      print_help(std::cerr, true);
      std::exit(0);
    }
    if (vars.count("help")) {
      print_help(std::cerr);
      std::exit(0);
    }

    po::notify(vars);
  } catch (const io_exception& ex) {
    std::cerr << ex.message() << std::endl << std::endl;
    print_help(std::cerr);
    std::exit(2);
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl << std::endl;
    print_help(std::cerr);
    std::exit(2);
  }

  return vars;
}

void Main::print_help(std::ostream& os, bool show_advanced) {
  os << boost::str(boost::format(m_usage) % m_name % get_version().make_string()) << std::endl;
  os << m_options << std::endl;
  if (show_advanced) {
    os << m_advanced_options << std::endl;
  }
}

pid_t Main::launch(const char* name, main_f fn, std::vector<const char*> args) {
  int retcode;
  auto pid = fork();
  switch (pid) {
    case -1:  // failed to fork
      SPLOG("Fork failed while launching: %s", name);
      exit(1);
      break;
    case 0:  // child
      log_change_name(name);
      retcode = main_child(name, fn, args);
      if (retcode != 0) {
        SPLOG("Failed to launch: %s", name);
        std::exit(1);
      }
      break;
    default:  // parent
      break;
  }
  return pid;
}

int Main::main_child(const char* name, main_f factory, std::vector<const char*> args) {
  prctl(PR_SET_PDEATHSIG, SIGTERM);

  try {
    setproctitle(name);

    auto module = factory();
    args.insert(args.begin(), name);
    auto vars = module->parse_args(args.size(), (char**)args.data());
    return module->run(vars);
  } catch (const io_exception& io) {
    SPLOG("Unhandled exception: %s\n", io.message().c_str());
    return 1;
  }
}

size_t Main::get_free_port() {
  ServerSocket svs(0);
  size_t port = svs.address().port();
  svs.close();

  SPLOG("Listening on port %lu", port);
  return port;
}

void Main::initialize_app(const std::string& ref_dir, const std::string& log_file) {
  // Ignore SIGINT so we can control the exit order
  signal(SIGINT, SIG_IGN);

  // Make temporary directory
  if (m_tmp_dir == "") {
    m_tmp_dir = "/tmp";
  }

  if (!fs::is_directory(m_tmp_dir)) {
    fs::create_directories(m_tmp_dir);
  }

  m_tmp_dir = m_tmp_dir + "/spiral_XXXXXX";
  char buf[m_tmp_dir.length() + 1];
  strcpy(buf, m_tmp_dir.c_str());
  char* r = mkdtemp(buf);
  if (r == NULL) {
    throw io_exception("Unable to make temp directory");
  }
  m_tmp_dir = fs::canonical(buf).native();

  path tdir(m_tmp_dir);
  tdir.mkdir();
  m_tmp_dir_made = true;

  if (not m_stats_file.empty()) {
    m_stats.save_to(m_stats_file);
    m_stats.add("date", pt::to_iso_extended_string(pt::microsec_clock::universal_time()) + "Z");
  }

#if GPERFTOOLS
  if (not m_cpuprofile_dir.empty()) {
    m_stats.save_cpuprofile_to(m_cpuprofile_dir);
  }
#endif

  // Set the appropriate number of threads (min 2, "auto" == 1 per cpu)
  set_thread_count(m_requested_threads);

  // Initiate logging
  // fprintf(stderr, "Setting log file %s/log.txt\n", m_tmp_dir.c_str());
#ifdef GPERFTOOLS
  log_init(getproctitle().c_str(), 2, true);
#else
  if (log_file == "") {
    m_logger_pid = setup_logger(m_tmp_dir + "/log.txt", m_debug_log);
    m_log_file = m_tmp_dir + "/log.txt";
  } else {
    m_logger_pid = setup_logger(log_file, m_debug_log);
    m_log_file = log_file;
  }
#endif

  SPLOG("%s", m_cmdline.c_str());
  SPLOG(" bg version: %s", get_version().make_string().c_str());
  SPLOG(" os release: %s", get_os_release().c_str());
  SPLOG("     kernel: %s", get_uname().c_str());
  SPLOG("       node: %s", get_nodename().c_str());
  SPLOG("        cpu: %d", std::thread::hardware_concurrency());
  SPLOG("    sys_mem: %lu GB", get_system_mem() / 1024 / 1024 / 1024);
  if (get_mem_limit() == std::numeric_limits<uint64_t>::max()) {
    SPLOG("  mem_limit: unlimited");
  } else {
    SPLOG("  mem_limit: %lu GB", get_mem_limit() / 1024 / 1024 / 1024);
  }
  SPLOG("   tmp_free: %lu GB on %s", fs::space(m_tmp_dir).available / 1024 / 1024 / 1024,
        m_tmp_dir.c_str());
  SPLOG("    threads: %lu", get_thread_count());

  // Set config
  Config::set("storage_root", tdir);
  Config::set("resources_root", tdir);
  Config::set("path_bulkdata", tdir);
  Config::set("temp_root", tdir);
  Config::set("reference_path", ref_dir);
  Config::set("task_timeout", 1200);
  Config::set("task_max_timeouts", 1);
  Config::set("task_update_interval", 2);
  std::vector<bind_info> baddrs;
  baddrs.resize(1);
  baddrs[0].port = get_free_port();
  baddrs[0].ip = "127.0.0.1";
  Config::set("taskdb_bind_list", json_wrap(baddrs));
  Config::set("taskdb_backup_period_in_seconds", 5);
}

void Main::launch_daemons() {
  // Launch worker procs
  m_taskdb_pid = launch("biograph_taskdb", taskdb_main);
  sleep(1);

  m_normal_pid =
      launch("biograph_manager", manager_main,
             {"--profile", "normal", "--num_procs", std::to_string(get_thread_count()).c_str()});
  m_himem_pid =
      launch("biograph_manager", manager_main, {"--profile", "himem", "--num_procs", "1"});
}

void Main::cleanup(bool success) {
  // Check if an override is present, so that cleanup() can be called at any time without checking.
  if (not needs_cleanup()) {
    return;
  }

  SPLOG("Shutting it down.");

  // Remove unimported reads on abort
  if (m_tmp_dir_made and not m_tmp_dir.empty() and fs::exists(m_tmp_dir)) {
    std::cerr << "Cleaning up...\n";
    fs::directory_iterator end_itr;
    for (fs::directory_iterator i(m_tmp_dir); i != end_itr; i++) {
      if (i->path().filename().string().find("_reads_") != std::string::npos) {
        fs::remove(i->path());
      }
    }
  }

  if (m_tmp_dir_made and success) {
    if (m_keep_tmp and fs::exists(m_tmp_dir)) {
      std::cerr << "Retaining temp directory " << m_tmp_dir << std::endl;
    } else {
      try {
        fs::remove_all(fs::path(m_tmp_dir));
      } catch (...) {
        // Ignore
      }
    }
    SPLOG("Finished");
  } else {
    if (!m_log_file.empty()) {
      SPLOG("There was a problem with this run. See %s for more details.", m_log_file.c_str());
      std::cerr << "There was a problem with this run. See " << m_log_file << " for more details.\n";
    }
    if (m_tmp_dir_made and not m_tmp_dir.empty()) {
      SPLOG("tmp-dir retained in %s", m_tmp_dir.c_str());
      std::cerr << "tmp-dir retained in " << m_tmp_dir << std::endl;
    }
  }
  // Moved pid closing to after log writing was done
  std::array<pid_t, 4> pids{{m_normal_pid, m_himem_pid, m_taskdb_pid, m_logger_pid}};
  for (size_t i = 0; i < pids.size(); i++) {
    if (pids[i]) {
      if (kill(pids[i], SIGTERM) == 0) {
        waitpid(pids[i], NULL, 0);
      }
    }
  }
}

int Main::main(const std::string& name, int argc, char* argv[]) {
  CHECK(spiral_initted()) << "Must call spiral_init from original main()";
  m_name = name;
  int ret = 1;

  try {
    auto vars = parse_args(argc, argv);
    ret = run(vars);
    cleanup(true);
  } catch (const io_exception& io) {
    std::cerr << "Error: " << io.message() << std::endl;
    cleanup(false);
  } catch (const std::runtime_error& er) {
    std::cerr << er.what() << std::endl;
    cleanup(true);
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << std::endl << std::endl;
    cleanup(false);
  }

  return ret;
}

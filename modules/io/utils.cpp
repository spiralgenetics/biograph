#include "modules/io/utils.h"
#include "base/base.h"
#include "modules/io/log.h"

#include <stdio.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/prctl.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <termios.h>
#include <unistd.h>
#include <array>
#include <cassert>
#include <cstdarg>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/uuid/detail/sha1.hpp>

#ifdef _WIN32
int vasprintf(char** sptr, const char* fmt, va_list argv) {
  int wanted = vsnprintf(*sptr = NULL, 0, fmt, argv);
  if ((wanted < 0) || ((*sptr = (char*)malloc(1 + wanted)) == NULL)) return -1;

  return vsprintf(*sptr, fmt, argv);
}
#endif

std::string printstring(const char* format, ...) {
  va_list vl;
  va_start(vl, format);

  char* buf;
  int len = vasprintf(&buf, format, vl);
  std::string r(buf, len);
  free(buf);
  va_end(vl);

  return r;
}

std::string time_to_RFC3339(const std::time_t time) {
  std::array<char, 100> buf;
  std::strftime(buf.data(), buf.size(), "%FT%TZ", std::localtime(&time));
  std::string now_str{buf.data()};
  CHECK(now_str.back() == 'Z');
  return now_str;
}

// Print a progress bar on stdout
void print_progress(float progress, int width) {
  // clamp to 100%
  if (progress > 1.0) {
    progress = 1.0;
  }
  if (isatty(fileno(stdout))) {
    std::cerr << "[";
    int pos = width * progress;
    for (int i = 0; i < width; ++i) {
      if (i < pos) {
        std::cerr << "=";
      } else if (i == pos) {
        std::cerr << ">";
      } else {
        std::cerr << " ";
      }
    }
    std::cerr << "] " << boost::format("%0.2f %%\r") % (progress * 100);
  } else {
    // Not a tty; output progress at most once every five minutes or 5% to avoid
    // filling up the log with progress updates.
    static time_t last_output = 0;
    static float prev_progress = 0.0;
    time_t now = std::time(0);
    if ((now >= (last_output + 300)) or (fabs(progress - prev_progress) > 0.05)) {
      std::cerr << "Progress = " << boost::format("%0.2f %%\n") % (progress * 100);
      last_output = now;
      prev_progress = progress;
    }
  }
  std::cerr.flush();
}

// Enable / disable echo on keystroke
void setecho(bool enable) {
  if (!isatty(fileno(stdin))) {
    return;
  }
  struct termios tty;
  tcgetattr(STDIN_FILENO, &tty);
  if (enable) {
    tty.c_lflag |= ECHO;
  } else {
    tty.c_lflag &= ~ECHO;
  }
  tcsetattr(STDIN_FILENO, TCSANOW, &tty);
}

// Return the current terminal width (min 80)
unsigned int get_terminal_width() {
  if (!isatty(fileno(stdin))) {
    return 80;
  }

  struct winsize win;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &win);

  if (win.ws_col < 80) {
    return 80;
  }
  return win.ws_col;
}

void set_mem_limit(uint64_t max_mem) {
  rlimit sys_limit;
  ::getrlimit(RLIMIT_AS, &sys_limit);

  if (max_mem == 0) {
    sys_limit.rlim_cur = sys_limit.rlim_max;
  } else {
    // std::cerr << "max_mem == " << max_mem << std::endl;
    sys_limit.rlim_cur = std::min(max_mem, sys_limit.rlim_max);
  }
  // std::cerr << "  set to: " << sys_limit.rlim_max << std::endl;

  // NOTE: RLIMIT_AS also counts against mmaps!
  if (::setrlimit(RLIMIT_AS, &sys_limit) != 0) {
    SPLOG("Unable to set memory limit: %s", strerror(errno));
  }
}

uint64_t get_mem_limit() {
  rlimit sys_limit;
  ::getrlimit(RLIMIT_AS, &sys_limit);

  return sys_limit.rlim_cur;
}

uint64_t get_system_mem() { return sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGE_SIZE); }

std::string get_uname() {
  struct utsname un;
  if (uname(&un) == -1) {
    return "uname() failed";
  }
  return boost::str(boost::format("%s %s %s %s") % un.sysname % un.release % un.version %
                    un.machine);
}

std::string get_nodename() {
  struct utsname un;
  if (uname(&un) == -1) {
    return "uname() failed";
  }
  return std::string(un.nodename);
}

std::string get_os_release() {
  if (!fs::exists("/etc/os-release")) {
    return "Unknown";
  }

  std::ifstream ifs("/etc/os-release");
  std::stringstream os_release;
  os_release << ifs.rdbuf();
  ifs.close();

  boost::regex re(".*^PRETTY_NAME=\"(.*?)\".*");
  boost::smatch match;
  std::string os_release_str = os_release.str();
  if (boost::regex_match(os_release_str, match, re)) {
    return match[1].str();
  }

  return "Unknown";
}

// Run a command through shell, return the output as a std::string
std::string easy_exec(const char* cmd) {
  char buffer[128];
  std::string result = "";
  std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (!feof(pipe.get())) {
    if (fgets(buffer, sizeof(buffer) / sizeof(buffer[0]), pipe.get()) != NULL) result += buffer;
  }
  return result;
}

std::string easy_exec(std::string cmd) { return easy_exec(cmd.c_str()); }

// tilde path expansion. Returns an empty string if HOME is not set.
std::string expand_home(const std::string& path) {
  auto new_path = std::string(path);
  if (path.substr(0, 2) != "~/") {
    return new_path;
  }

  const char* home = std::getenv("HOME");
  if (home == nullptr) {
    return "";
  }

  return new_path.replace(0, 1, std::string(home));
}

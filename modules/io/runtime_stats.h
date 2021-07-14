#pragma once

#include "json_spirit.h"
#include "modules/io/version.h"

#include <sys/resource.h>
#include <sys/time.h>
#include <boost/filesystem.hpp>
#include <chrono>

namespace js = json_spirit;

struct runtime_stage {
  std::string name;
  time_t end_time = 0;
  double wall_seconds = 0;
  double cpu_user_seconds = 0;
  double cpu_sys_seconds = 0;
};

typedef std::vector<runtime_stage> stages_t;

class runtime_stats {
 public:
  // default constructor
  runtime_stats();

  // construct with a file path
  runtime_stats(const std::string& out_file);

  // save on destruct
  ~runtime_stats() { save(); }

  // Add a stat
  void add(const js::Pair& stat);

  template <typename T>
  void add(const std::string& name, T const& value) {
    add(js::Pair(name, value));
  };

  // Add a timing
  runtime_stage add_timing(const std::string& name, const std::time_t& time);
  runtime_stage add_timing(const std::string& name) { return add_timing(name, std::time(0)); }

  // Getters
  const stages_t& get_stages() const { return m_stages; }
  const js::Object& get_stats() const { return m_stats; }

  // Clear all stats and timings
  void clear() {
    m_stats.clear();
    m_stages.erase(m_stages.begin()++);
  }

  // Write the stats to m_out_file. Does nothing if m_out_file is not set.
  void save();

  // Change the location of m_out_file. Does not actually save.
  void save_to(const std::string& out_file) { m_out_file = out_file; }

  // Log stage start, stage end, and save timings.
  void start_stage(const std::string& stage_name);
  void end_stage(const std::string& stage_name);

#if GPERFTOOLS
  // Change the location of the cpu profile directory to be used for subsequent stages.
  void save_cpuprofile_to(const std::string& cpuprofile_dir);
#endif

 private:
  std::string m_out_file;
  js::Object m_stats;
  stages_t m_stages;
  std::string m_cur_stage;
  std::chrono::time_point<std::chrono::steady_clock> m_stage_start_time;
  struct rusage m_stage_start_rusage;
#if GPERFTOOLS
  std::string m_cpuprofile_dir;
#endif
};

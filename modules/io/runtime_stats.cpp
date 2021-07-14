#include "modules/io/runtime_stats.h"
#include "base/base.h"
#include "modules/io/io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"

#include <fstream>

#if GPERFTOOLS
#include <gperftools/profiler.h>
#endif

namespace js = json_spirit;

const std::string k_start_timings = "runtime_stats_start";
const std::string k_end_timings = "runtime_stats_end";

// default constructor
runtime_stats::runtime_stats() {
  add_timing(k_start_timings);
};

// construct with a file path
runtime_stats::runtime_stats(const std::string& out_file) {
  m_out_file = out_file;
  add_timing(k_start_timings);
}

void runtime_stats::add(const js::Pair& stat) {
  m_stats.push_back(stat);
}

namespace {

double timeval_diff(struct timeval start, struct timeval end) {
  double secs = end.tv_sec - start.tv_sec;
  secs += (int64_t(end.tv_usec) - int64_t(start.tv_usec)) / 1000000.;
  return secs;
}

}  // namespace

runtime_stage runtime_stats::add_timing(const std::string& name, const std::time_t& now) {
#if GPERFTOOLS
   if (not m_cpuprofile_dir.empty()) {
     ProfilerStop();
     boost::filesystem::rename(m_cpuprofile_dir + "/final.prof",
                               m_cpuprofile_dir + "/" + name + ".prof");
     ProfilerStart((m_cpuprofile_dir + "/final.prof").c_str());
   }
#endif
  runtime_stage t;
  t.name = name;
  t.end_time = now;
  if (m_cur_stage == name) {
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = end_time - m_stage_start_time;
    t.wall_seconds = time_diff.count();

    struct rusage stage_end_rusage;
    if (getrusage(RUSAGE_SELF, &stage_end_rusage) < 0) {
      throw(io_exception("getrusage failed"));
    }

    t.cpu_user_seconds = timeval_diff(m_stage_start_rusage.ru_utime, stage_end_rusage.ru_utime);
    t.cpu_sys_seconds = timeval_diff(m_stage_start_rusage.ru_stime, stage_end_rusage.ru_stime);
  }
  m_stages.push_back(t);
  return t;
}

void runtime_stats::save() {
  if (m_out_file.empty()) {
    return;
  }

  // The end is the last time we save.
  if (m_stages.back().name == k_end_timings) {
    m_stages.erase(m_stages.end());
  }
  add_timing(k_end_timings, std::time(0));

  js::Object stats = m_stats;
  stats.push_back(js::Pair("runtime_seconds", m_stages.back().end_time - m_stages.front().end_time));

  js::Array json_stages;
  for (const auto& stage : get_stages()) {
    if (stage.name == k_start_timings || stage.name == k_end_timings) {
      continue;
    }
    js::Object json_stage;
    json_stage.push_back(js::Pair("name", stage.name));
    if (stage.wall_seconds || stage.cpu_user_seconds || stage.cpu_sys_seconds) {
      json_stage.push_back(js::Pair("wall_seconds", stage.wall_seconds));
      json_stage.push_back(js::Pair("cpu_user_seconds", stage.cpu_user_seconds));
      json_stage.push_back(js::Pair("cpu_sys_seconds", stage.cpu_sys_seconds));
    }
    json_stages.push_back(json_stage);
  }
  stats.push_back(js::Pair("stages", json_stages));

  std::ofstream os(m_out_file);
  js::write(stats, os);
  if (not os.good()) {
    throw io_exception("Could not write stats to " + m_out_file);
  }
  os.close();
}

void runtime_stats::start_stage(const std::string& stage_name) {
  CHECK_EQ(m_cur_stage, "") << "Cannot start stage " << stage_name << " with " << m_cur_stage
                            << " already in progress";
  if (getrusage(RUSAGE_SELF, &m_stage_start_rusage) < 0) {
    throw(io_exception("getrusage falied"));
  }
  m_stage_start_time = std::chrono::steady_clock::now();
  m_cur_stage = stage_name;
  SPLOG("Start Stage::%s", stage_name.c_str());
}

void runtime_stats::end_stage(const std::string& stage_name) {
  if (m_cur_stage.empty()) {
    CHECK(!m_cur_stage.empty()) << "Expecting to end stage " << stage_name;
  } else {
    CHECK_EQ(m_cur_stage, stage_name) << "Expecting to end stage " << stage_name
                                      << " but we are in stage " << m_cur_stage;
  }
  auto t = add_timing(m_cur_stage);
  m_cur_stage.clear();
  double tot_cpu = t.cpu_user_seconds + t.cpu_sys_seconds;
  double avg_parallel = 0;

  if (t.wall_seconds) {
    avg_parallel = tot_cpu / t.wall_seconds;
  }
  std::string human_wall_seconds;
  if (t.wall_seconds >= 3600) {
    human_wall_seconds = printstring(" (%dh%02dm%02ds)", int(t.wall_seconds) / 3600,
                                     (int(t.wall_seconds) % 3600) / 60, int(t.wall_seconds) % 60);
  } else if (t.wall_seconds > 60) {
    human_wall_seconds =
        printstring(" (%dm%ds)", int(t.wall_seconds) / 60, int(t.wall_seconds) % 60);
  }
  double user_percent = 0;
  if (tot_cpu) {
    user_percent = t.cpu_user_seconds * 100 / tot_cpu;
  }
  SPLOG(
      "End Stage::%s  Wall time: %.2f sec%s Avg parallelism: %.2f  CPU time: %.2f sec:  %.2fs sec "
      "user (%.2f%%) + %.2f sec system)",
      stage_name.c_str(), t.wall_seconds, human_wall_seconds.c_str(), avg_parallel, tot_cpu,
      t.cpu_user_seconds, user_percent, t.cpu_sys_seconds);
}

#if GPERFTOOLS
void runtime_stats::save_cpuprofile_to(const std::string& cpuprofile_dir) {
  m_cpuprofile_dir = cpuprofile_dir;
  if (not boost::filesystem::exists(cpuprofile_dir)) {
    boost::filesystem::create_directory(cpuprofile_dir);
  }

  ProfilerStart((cpuprofile_dir + "/final.prof").c_str());
}
#endif

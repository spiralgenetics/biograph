#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/reversable_tracer.h"
#include "modules/variants/scaffold.h"

namespace variants {

class trace_ref {
 public:
  trace_ref(const assemble_options& options, scaffold_pipeline_interface* output_f);
  ~trace_ref();

  void add_entire_reference();
  void add_scaffold(const std::string& scaffold_name);
  void add_scaffold_range(const std::string& scaffold_name, size_t start, size_t limit);

  assemble_stats assemble(progress_handler_t progress = null_progress_handler);

  // True if no work has been queued.
  bool empty() const;

  // Abort the current trace in progress.
  void abort_trace();

  static scaffold ref_to_scaffold(const reference* ref, const std::string& scaffold_name);

  static void display_in_progress();

  static bool g_verbose_trace_work;

 private:
  struct work_info {
    ~work_info();

    bool skip_fwd;
    bool skip_rev;
    aoffset_t start;
    aoffset_t limit;
    std::unique_ptr<pipeline_interface> p;
    std::unique_ptr<reversable_tracer> pop;
    std::unique_ptr<reversable_tracer> rc_pop;
    std::string scaffold_name;
    std::shared_ptr<scaffold> s;

    std::function<void(const assembly&, bool /* anchored on right */)> report_anchor_drop_func;

    std::string to_string() const {
      return printstring("%s[%d,%d)", scaffold_name.c_str(), start, limit);
    }
  };

  static void note_work_start(const std::unique_ptr<work_info>& w, const std::string& work_desc);
  static void note_work_finish(const std::unique_ptr<work_info>& w, const std::string& work_desc);
  static std::string work_in_progress();
  void abort_work(std::unique_ptr<work_info> w);

  static std::mutex g_in_progress_mu;
  using in_progress_key_t = std::pair<work_info*, std::string /* desc */>;
  static std::map<in_progress_key_t, time_t /* start time */> g_in_progress;
  std::shared_ptr<scaffold> get_scaffold(const std::string& scaffold_name) const;
  assemble_stats execute_work(std::unique_ptr<work_info> w) const;
  assemble_stats execute_work_direction(work_info* w, bool rev_comp,
                                        const assemble_options& opts) const;

  assemble_options m_options;
  std::vector<std::unique_ptr<work_info>> m_work;
  scaffold_pipeline_interface* m_output_f = nullptr;

  bool m_aborted = false;
};

}  // namespace variants

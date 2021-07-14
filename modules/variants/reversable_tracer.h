#pragma once

#include "modules/variants/pair_counter.h"
#include "modules/variants/tracer.h"
#include "modules/variants/pop_tracer.h"
#include "modules/variants/discovery/state.h"

namespace variants {

class reversable_tracer {
 public:
  reversable_tracer(bool rev_comp, const assemble_options& options);

  assemble_stats assemble(assemble_pipeline_interface* output,
                          progress_handler_t progress = null_progress_handler);
  assemble_stats assemble(aoffset_t start_offset, aoffset_t limit_offset,
                          assemble_pipeline_interface* output,
                          progress_handler_t = null_progress_handler);
  void output_path_debug_dot(std::ostream& os) const {
    m_tracer->output_path_debug_dot(os);
  }

  // If in pop tracer mode, wraps the given half aligned assembly reporter to
  // gather data the pop tracer needs.
  std::function<void(const assembly&, bool /* right anchor */)>
  wrap_report_anchor_drop_for_pop_tracer(
      const std::function<void(const assembly&, bool /* right anchor */)>& in);

  // Adds a potential read in an approximate location.  If "rev_comp" is true, this read is already
  // reversed.
  void add_approx_read(uint32_t read_id, aoffset_t start_limit, aoffset_t end_limit, bool rev_comp);

 private:
  bool m_rev_comp;
  aoffset_t m_ref_end_pos;

  scaffold m_rev_scaffold;
  assemble_options m_options;
  boost::optional<discovery::state> m_bidir_tracer;
  boost::optional<tracer> m_tracer;
  boost::optional<pop_tracer> m_pop_tracer;

  std::unique_ptr<pair_counter> m_pair_counter;
};

}  // namespace variants

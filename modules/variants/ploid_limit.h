#pragma once

#include "modules/variants/assemble.h"

namespace variants {

// Filters assemblies to the number of allowed alleles based on score.
//
// TODO(nils): Rename and redocument, since now this more merges
// matching assemblies than it limits ploids.
class ploid_limiter : public sorted_output_pipeline_step {
 public:
  ploid_limiter(const assemble_options& options, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output)), m_options(options) {
    set_expected_order(assembly::left_offset_less_than);
  }
  ploid_limiter() = delete;
  ~ploid_limiter();

  void on_assembly(assembly_ptr a) override;

  // TODO(nils): Ploid limiter shoold just merge, not limit, maybe?
  // Figure out if we actually want to do this and update tests to not
  // test for limiting.
  void set_max_ploids(unsigned max_ploids) { m_max_ploids = max_ploids; }

 private:
  acost_t calc_score(const assembly& a) const;
  void output_active();
  void flush_queued();
  void do_deploid();
  void ploid_flush();

  assemble_options m_options;
  pipeline_step_t m_output;

  aoffset_t m_cur_offset = 0;

  // Limit output to this many ploids.  Ideally this would be
  // unlimited so the genotyper can do the full coverage calculation,
  // but we should limit it some here for performance reasons.g0
  unsigned m_max_ploids = 20;

  // Number of non-reference assemblies present in m_active.
  unsigned m_var_active = 0;
  std::multimap<aoffset_t /* right offset */, assembly_ptr> m_active;

  // When n_active > ploid_limit, we queue to m_deploid_scores.  When
  // n_active <= ploid_limit again, we execute the sort-and-deploid
  // stage.
  std::multimap<acost_t /* score */, assembly_ptr> m_deploid_scores;
};

}  // namespace variants

#pragma once

#include "modules/variants/assemble.h"

namespace variants {

// Runs all the standard pieces of a variant calling pipeline.
// This could be used, for instance, in between a reversable tracer and a
// VCF exporter.
class assemble_pipeline : public pipeline_interface {
 public:
  assemble_pipeline(const assemble_options& options, pipeline_step_t output);
  assemble_pipeline() = delete;
  ~assemble_pipeline();

  pipeline_step_t make_parallel_input() override;

  template <typename T, typename... Args>
  void add_step(Args... args) {
    add_step_internal<T, Args...>(m_ser, args...);
  }
  void add_standard_variants_pipeline();

 private:
  template <typename T, typename... Args>
  void add_step_internal(std::vector<std::function<pipeline_step_t(pipeline_step_t)>>& step_makers,
                         Args... args) {
    std::function<pipeline_step_t(Args..., pipeline_step_t)> step_constructor =
        make_unique<T, Args..., pipeline_step_t>;
    step_makers.push_back(std::bind(step_constructor, args..., std::placeholders::_1));
  }

  pipeline_step_t make_pipeline(
      const std::vector<std::function<pipeline_step_t(pipeline_step_t)>>& steps,
      pipeline_step_t output);

 private:
  std::mutex m_mu;

  assemble_options m_options;
  pipeline_step_t m_first_ser_step;
  pipeline_step_t m_output;
  // Parallel steps; up to and including the sort.
  std::vector<std::function<pipeline_step_t(pipeline_step_t)>> m_par;
  // Serial steps; after the sort.
  std::vector<std::function<pipeline_step_t(pipeline_step_t)>> m_ser;
};

}  // namespace variants

#pragma once

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "modules/variants/assemble.h"

namespace variants {

class rvg_exclude : public sorted_output_pipeline_step {
 public:
  rvg_exclude(const assemble_options& opts, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output)), m_options(opts) {}
  rvg_exclude() = delete;
  ~rvg_exclude() { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  void flush();

  assemble_options m_options;

  absl::flat_hash_map<size_t /* assembly id */, std::vector<assembly_ptr>> m_backlog;
  absl::flat_hash_set<size_t /* assembly id */> m_known_inphase;
};

}  // namespace variants

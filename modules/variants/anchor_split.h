#pragma once

#include "modules/variants/assemble.h"

namespace variants {

// Splits the ref-matching anchors off all assemblies, so each
// assembly becomes 2 reference-only assemblies and a
// non-reference-only assembly.  If the "trace_reference_assemblies"
// option isn't set, the reference-only assemblies are skipped.
class anchor_splitter : public sorted_output_pipeline_step {
 public:
  anchor_splitter(const assemble_options& options, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output)), m_options(options) {
    set_expected_order(assembly::left_offset_less_than);
  }
  anchor_splitter() = delete;

  void on_assembly(assembly_ptr a) override;

 private:
  assemble_options m_options;
};

}  // namespace variants

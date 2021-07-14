#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

// Trims all sections of incoming assemblies that match reference.

class ref_trimmer : public sorted_output_pipeline_step {
 public:
  ref_trimmer(const assemble_options& options, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output), true /* old sort order */),
        m_options(options),
        m_scaffold(options.scaffold) {
    CHECK(m_scaffold);
  }
  ref_trimmer() = delete;
  ~ref_trimmer() override { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  void flush();

  assemble_options m_options;
  const scaffold* m_scaffold = nullptr;

  // Maximum number of bases to backtrack when trimming reference
  // bases on the right side of a half-anchored assembly.
  static constexpr aoffset_t k_max_backtrack_len = 1000;
};

}  // namespace variants

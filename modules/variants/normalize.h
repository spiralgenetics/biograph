#pragma once

#include "modules/variants/assemble.h"

namespace variants {

class normalizer : public sorted_output_pipeline_step {
 public:
  normalizer(const assemble_options& options, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output), true /* old sort order */), m_options(options) {
    CHECK(m_options.scaffold);
  }
  normalizer() = delete;
  ~normalizer() { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  void advance_to(aoffset_t offset);
  void flush();

  assemble_options m_options;
  aoffset_t m_cur_offset = 0;
};

// Pad things that have empty sequence or reference by adding the next
// reference base to the left.  This is how VCF expects things.
class vcf_padder : public sorted_output_pipeline_step {
 public:
  vcf_padder(const assemble_options& options, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output), true /* old sort order */), m_options(options) {
    CHECK(m_options.scaffold);
  }
  vcf_padder() = delete;
  ~vcf_padder() { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  void flush();

  assemble_options m_options;
};

}  // namespace variants

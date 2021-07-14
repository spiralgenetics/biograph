#pragma once

#include "modules/variants/assemble.h"

namespace variants {

// Genotypes input assemblies based on the min depth of the variants.
class simple_genotype_filter : public assemble_pipeline_interface {
 public:
  simple_genotype_filter(const assemble_options& options, pipeline_step_t output)
      : m_options(options), m_output(std::move(output)) {
    CHECK(options.scaffold);
  }
  ~simple_genotype_filter() override = default;

  void on_assembly(assembly_ptr a) override;

 private:
  void discard(assembly_ptr a);

  assemble_options m_options;
  pipeline_step_t m_output;
};

}  // namespace variants

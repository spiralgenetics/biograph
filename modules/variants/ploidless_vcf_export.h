#pragma once

#include "modules/bio_base/reference.h"
#include "modules/io/io.h"
#include "modules/variants/assemble.h"

namespace variants {

class ploidless_vcf_export : public assemble_pipeline_interface {
 public:
  ploidless_vcf_export(const assemble_options& options, std::string scaffold_name,
                       const std::function<void(const std::string&)>& output_f)
      : m_options(options), m_scaffold_name(scaffold_name), m_output(output_f) {}

  void on_assembly(assembly_ptr a) override;

  static std::string header(const assemble_options& options,
                            const std::map<std::string, std::string>& extra_headers,
                            const std::string sample_name = "SAMPLE");

 private:
  assemble_options m_options;
  std::string m_scaffold_name;
  std::function<void(const std::string&)> m_output;
};

}  // namespace variants

#include "python/biograph/variants/dedup_cov_reads.h"

#include "modules/variants/dedup_cov_reads.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) dedup_cov_reads_generator
    : public asm_pipeline_wrapper {
 public:
  dedup_cov_reads_generator(const reference_wrapper& ref, const std::string& scaffold_name,
                            object input);

  void init();
  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  assemble_options m_options;
  ::variants::scaffold m_scaffold;

  std::unique_ptr<dedup_cov_reads> m_dedup_cov_reads;
};

dedup_cov_reads_generator::dedup_cov_reads_generator(const reference_wrapper& ref,
                                                     const std::string& scaffold_name, object input)
    : asm_pipeline_wrapper(input) {
  m_options.scaffold_name = scaffold_name;
  m_scaffold = trace_ref::ref_to_scaffold(ref.get_reference().get(), scaffold_name);
  m_options.scaffold = &m_scaffold;
}

void dedup_cov_reads_generator::init() {
  m_dedup_cov_reads.reset(new dedup_cov_reads(m_options, make_pipeline_output()));
}

void dedup_cov_reads_generator::on_input(assembly_ptr a) { m_dedup_cov_reads->add(std::move(a)); }

void dedup_cov_reads_generator::on_input_done() { m_dedup_cov_reads.reset(); }

void bind_dedup_cov_reads(module& m) {
  define_pipeline_generator<dedup_cov_reads_generator, const reference_wrapper& /* ref */,
                            const std::string& /* scaffold_name */, object /* input*/>(
      m, "dedup_cov_reads", "DedupCovReadsGenerator", arg("ref"), arg("scaffold_name"),
      arg("input"), R"DOC(how does this work?)DOC");
}

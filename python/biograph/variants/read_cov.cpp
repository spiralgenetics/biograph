#include "python/biograph/variants/read_cov.h"

#include "modules/variants/add_ref.h"
#include "modules/variants/read_cov.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) read_cov_generator : public par_asm_pipeline_wrapper {
 public:
  read_cov_generator(const std::shared_ptr<readmap>& rm, object input, int max_reads_per_entry,
                     int max_coverage_paths);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  assemble_options m_options;

  std::shared_ptr<readmap> m_readmap;

  std::unique_ptr<read_cov> m_read_cov;
};

read_cov_generator::read_cov_generator(const std::shared_ptr<readmap>& rm, object input,
                                       int max_reads_per_entry, int max_coverage_paths)
    : par_asm_pipeline_wrapper(input), m_readmap(rm) {
  m_options.read_cov_max_reads_per_entry = max_reads_per_entry;
  m_options.max_coverage_paths = max_coverage_paths;
  m_options.seqset = rm->get_seqset().get();
  m_options.readmap = rm.get();
}

void read_cov_generator::init() {
  m_read_cov.reset(new read_cov(m_options, make_pipeline_output()));
}

void read_cov_generator::on_input(assembly_ptr a) { m_read_cov->add(std::move(a)); }

void read_cov_generator::on_input_done() { m_read_cov.reset(); }

void bind_read_cov(module& m) {
  define_pipeline_generator<read_cov_generator, const std::shared_ptr<readmap>&, object /* input */,
                            int /* max_reads_per_entry */, int /* max_coverage_paths */>(
      m, "generate_read_cov", "ReadCovGenerator", arg("rm"), arg("input"),
      arg("max_reads_per_entry") = 0, arg("max_coverage_paths") = 0,
      R"DOC(how does this work?)DOC");
}

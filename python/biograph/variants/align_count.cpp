#include "python/biograph/variants/align_count.h"

#include "modules/variants/add_ref.h"
#include "modules/variants/align_count.h"
#include "modules/variants/pair_cov.h"
#include "modules/variants/read_cov.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) align_count_generator
    : public par_asm_pipeline_wrapper {
 public:
  align_count_generator(object input);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  assemble_options m_options;
  ::variants::scaffold m_scaffold;

  std::unique_ptr<align_count> m_cov;
  std::shared_ptr<readmap> m_readmap;
};

align_count_generator::align_count_generator(object input) : par_asm_pipeline_wrapper(input) {}

void align_count_generator::init() {
  m_cov.reset(new align_count(m_options, make_pipeline_output()));
}

void align_count_generator::on_input(assembly_ptr a) { m_cov->add(std::move(a)); }

void align_count_generator::on_input_done() { m_cov.reset(); }

void bind_align_count(module& m) {
  define_pipeline_generator<align_count_generator, object /* input */>(
      m, "align_count", "AlignCountGenerator", (arg("input")), R"DOC(how does this work?)DOC");
}

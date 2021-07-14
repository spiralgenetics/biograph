#include "python/biograph/variants/pair_edge_cov.h"

#include "modules/variants/add_ref.h"
#include "modules/variants/pair_cov.h"
#include "modules/variants/pair_edge_cov.h"
#include "modules/variants/read_cov.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) pair_edge_cov_generator
    : public par_asm_pipeline_wrapper {
 public:
  pair_edge_cov_generator(object input);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  assemble_options m_options;
  ::variants::scaffold m_scaffold;

  std::unique_ptr<pair_edge_cov> m_cov;
  std::shared_ptr<readmap> m_readmap;
};

pair_edge_cov_generator::pair_edge_cov_generator(object input) : par_asm_pipeline_wrapper(input) {}

void pair_edge_cov_generator::init() {
  m_cov.reset(new pair_edge_cov(m_options, make_pipeline_output()));
}

void pair_edge_cov_generator::on_input(assembly_ptr a) { m_cov->add(std::move(a)); }

void pair_edge_cov_generator::on_input_done() { m_cov.reset(); }

void bind_pair_edge_cov(module& m) {
  define_pipeline_generator<pair_edge_cov_generator, object /* input */>(
      m, "generate_pair_edge_cov", "PairEdgeCovGenerator", (arg("input")),
      R"DOC(how does this work?)DOC");
}

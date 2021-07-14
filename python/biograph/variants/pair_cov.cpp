#include "python/biograph/variants/pair_cov.h"

#include "modules/variants/add_ref.h"
#include "modules/variants/pair_cov.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) pair_cov_generator : public par_asm_pipeline_wrapper {
 public:
  pair_cov_generator(const std::shared_ptr<readmap>& rm, object input,
                     const unsigned min_pair_distance, const unsigned max_pair_distance);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  assemble_options m_options;

  std::shared_ptr<readmap> m_readmap;
  std::unique_ptr<pair_cov> m_pair_cov;
};

pair_cov_generator::pair_cov_generator(const std::shared_ptr<readmap>& rm, object input,
                                       const unsigned min_pair_distance,
                                       const unsigned max_pair_distance)
    : par_asm_pipeline_wrapper(input), m_readmap(rm) {
  m_options.seqset = rm->get_seqset().get();
  m_options.readmap = rm.get();

  m_options.min_pair_distance = min_pair_distance;
  m_options.max_pair_distance = max_pair_distance;
}

void pair_cov_generator::init() {
  m_pair_cov.reset(new pair_cov(m_options, make_pipeline_output()));
}

void pair_cov_generator::on_input(assembly_ptr a) { m_pair_cov->add(std::move(a)); }

void pair_cov_generator::on_input_done() { m_pair_cov.reset(); }

void bind_pair_cov(module& m) {
  define_pipeline_generator<pair_cov_generator, const std::shared_ptr<readmap>&, object /* input */,
                            unsigned /* min_pair_distance */, unsigned /* max_pair_distance */>(
      m, "generate_pair_cov", "PairCovGenerator", arg("rm"), arg("input"),
      arg("min_insert_size") = 200, arg("max_insert_size") = 1000, R"DOC(how does this work?)DOC");
}

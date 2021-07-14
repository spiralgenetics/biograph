#include "python/biograph/variants/place_pair_cov.h"

#include "absl/strings/match.h"
#include "absl/strings/str_split.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/place_pair_cov.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

#include <pybind11/functional.h>

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) place_pair_cov_generator
    : public par_asm_pipeline_wrapper {
 public:
  place_pair_cov_generator(object input, std::shared_ptr<readmap> rm, int min_insert_size,
                           int max_insert_size, int ideal_insert_size, int max_ambig);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  std::unique_ptr<place_pair_cov> m_place_pair_cov;
  assemble_options m_opts;
  place_pair_options m_popts;
};

place_pair_cov_generator::place_pair_cov_generator(object input, std::shared_ptr<readmap> rm,
                                                   int min_insert_size, int max_insert_size,
                                                   int ideal_insert_size, int max_ambig)
    : par_asm_pipeline_wrapper(input) {
  m_opts.seqset = rm->get_seqset().get();
  m_opts.readmap = rm.get();
  m_opts.min_pair_distance = min_insert_size;
  m_opts.max_pair_distance = max_insert_size;

  m_popts.max_ambig = max_ambig;
  m_popts.ideal_pair_distance = ideal_insert_size;
}

void place_pair_cov_generator::init() {
  m_place_pair_cov.reset(new place_pair_cov(m_opts, m_popts, make_pipeline_output()));
}

void place_pair_cov_generator::on_input(assembly_ptr a) { m_place_pair_cov->add(std::move(a)); }

void place_pair_cov_generator::on_input_done() {
  m_place_pair_cov->flush();
  m_place_pair_cov.reset();
}

void bind_place_pair_cov(module& m) {
  define_pipeline_generator<place_pair_cov_generator, object /* input */, std::shared_ptr<readmap>,
                            int, int, int, int>(
      m, "place_pair_cov", "PlacePairCovGenerator", arg("input"), arg("rm"),
      arg("min_insert_size") = 200, arg("max_insert_size") = 1000, arg("ideal_insert_size") = 400,
      arg("max_ambig") = 15, R"DOC(how does this work?)DOC");
}

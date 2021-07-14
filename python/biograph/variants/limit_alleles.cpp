#include "python/biograph/variants/limit_alleles.h"

#include "modules/variants/limit_alleles.h"
#include "modules/variants/phase.h"
#include "modules/variants/scaffold.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/pipeline.h"
#include "python/common.h"

#include <pybind11/functional.h>

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) limit_alleles_generator : public asm_pipeline_wrapper {
 public:
  limit_alleles_generator(size_t max_alleles, const limit_alleles::sort_func_t& sort_func,
                          const limit_alleles::on_limited_func_t& on_limited, object input);

  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  unsigned m_max_alleles;
  limit_alleles::sort_func_t m_sort_func;
  limit_alleles::on_limited_func_t m_on_limited_func;

  std::unique_ptr<limit_alleles> m_limit_alleles;
};

limit_alleles_generator::limit_alleles_generator(size_t max_alleles,
                                                 const limit_alleles::sort_func_t& sort_func,
                                                 const limit_alleles::on_limited_func_t& on_limited,
                                                 object input)
    : asm_pipeline_wrapper(input),
      m_max_alleles(max_alleles),
      m_sort_func(sort_func),
      m_on_limited_func(on_limited) {}

void limit_alleles_generator::init() {
  m_limit_alleles.reset(
      new limit_alleles(m_max_alleles, m_sort_func, m_on_limited_func, make_pipeline_output()));
}

void limit_alleles_generator::on_input(assembly_ptr a) { m_limit_alleles->add(std::move(a)); }

void limit_alleles_generator::on_input_done() { m_limit_alleles.reset(); }

void bind_limit_alleles(module& m) {
  define_pipeline_generator<limit_alleles_generator, size_t /* max_alleles */,
                            limit_alleles::sort_func_t, limit_alleles::on_limited_func_t,
                            object /* input */>(
      m, "limit_alleles", "LimitAllelesGenerator", arg("max_alleles"), arg("sort_alleles"),
      arg("on_limited"), arg("input"), R"DOC(how does this work?)DOC");
}

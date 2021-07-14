#include "python/biograph/variants/filter_dup_align.h"

#include "modules/variants/add_ref.h"
#include "modules/variants/filter_dup_align.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

#include <pybind11/functional.h>

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) filter_dup_align_generator
    : public par_asm_pipeline_wrapper {
 public:
  filter_dup_align_generator(const filter_dup_align::sort_func_t& sort_func, object input);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  filter_dup_align::sort_func_t m_sort_func;
  std::unique_ptr<filter_dup_align> m_filter_dup_align;
};

filter_dup_align_generator::filter_dup_align_generator(
    const filter_dup_align::sort_func_t& sort_func, object input)
    : par_asm_pipeline_wrapper(input), m_sort_func(sort_func) {}

void filter_dup_align_generator::init() {
  m_filter_dup_align.reset(new filter_dup_align(
      [this](std::vector<assembly_ptr> input) -> std::vector<assembly_ptr> {
        gil_scoped_acquire gil;
        return m_sort_func(std::move(input));
      },
      make_pipeline_output()));
}

void filter_dup_align_generator::on_input(assembly_ptr a) { m_filter_dup_align->add(std::move(a)); }

void filter_dup_align_generator::on_input_done() { m_filter_dup_align.reset(); }

void bind_filter_dup_align(module& m) {
  define_pipeline_generator<filter_dup_align_generator, filter_dup_align::sort_func_t,
                            object /* input */>(m, "filter_dup_align", "FilterDupAlignGenerator",
                                                arg("sort_func"), arg("input"),
                                                R"DOC(how does this work?)DOC");
}

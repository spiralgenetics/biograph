#include "python/biograph/variants/align_reads.h"

#include "modules/variants/align_reads.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

#include <pybind11/functional.h>

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) align_reads_generator
    : public par_asm_pipeline_wrapper {
 public:
  align_reads_generator(object input, const align_reads::on_aligned_func_t& f, bool refskip_anchor);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  align_reads::on_aligned_func_t m_on_aligned;

  std::unique_ptr<align_reads> m_align;
  bool m_refskip_anchor = false;
};

align_reads_generator::align_reads_generator(object input,
                                             const align_reads::on_aligned_func_t& on_aligned,
                                             bool refskip_anchor)
    : par_asm_pipeline_wrapper(input), m_on_aligned(on_aligned), m_refskip_anchor(refskip_anchor) {}

void align_reads_generator::init() {
  m_align.reset(new align_reads(
      [this](read_id_set read_ids, aligned_read read) {
        gil_scoped_acquire gil;
        m_on_aligned(std::move(read_ids), std::move(read));
      },
      m_refskip_anchor, make_pipeline_output()));
}

void align_reads_generator::on_input(assembly_ptr a) { m_align->add(std::move(a)); }

void align_reads_generator::on_input_done() { m_align.reset(); }

void bind_align_reads(module& m) {
  class_<aligned_read>(m, "AlignedRead")
      .def_readwrite("left_offset", &aligned_read::left_offset)
      .def_readwrite("right_offset", &aligned_read::right_offset)
      .def_readwrite("cigar", &aligned_read::cigar)
      .def_readwrite("seq", &aligned_read::seq)
      .def("__str__", str_from_ostream<aligned_read>);

  define_pipeline_generator<align_reads_generator, object /* input */,
                            align_reads::on_aligned_func_t, bool /* refskip anchor */>(
      m, "align_reads", "AlignReadsGenerator", arg("input"), arg("on_aligned"),
      arg("refskip_anchor") = false, R"DOC(how does this work?)DOC");
}

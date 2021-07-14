#include "python/biograph/variants/apply_edges.h"

#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) apply_edges_py : public apply_edges_step {
 public:
  apply_edges_py(pipeline_step_t output, object on_assembly_edges)
      : apply_edges_step(std::move(output)), m_on_assembly_edges(on_assembly_edges) {}
  ~apply_edges_py() { flush(); };

  void on_assembly_edges(optional_aoffset reference_pos,
                         const std::vector<assembly_ptr>& left_edges,
                         const std::vector<assembly_ptr>& inserts,
                         const std::vector<assembly_ptr>& right_edges) override {
    gil_scoped_acquire gil;

    m_on_assembly_edges(reference_pos, left_edges, inserts, right_edges);
  }

 private:
  object m_on_assembly_edges;
};

class __attribute__((visibility("hidden"))) apply_edges_generator
    : public par_asm_pipeline_wrapper {
 public:
  apply_edges_generator(object input, object on_assembly_edges);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  object m_on_assembly_edges;
  std::unique_ptr<apply_edges_py> m_apply_edges;
};

apply_edges_generator::apply_edges_generator(object input, object on_assembly_edges)
    : par_asm_pipeline_wrapper(input), m_on_assembly_edges(on_assembly_edges) {}

void apply_edges_generator::init() {
  m_apply_edges.reset(new apply_edges_py(make_pipeline_output(), m_on_assembly_edges));
}

void apply_edges_generator::on_input(assembly_ptr a) { m_apply_edges->add(std::move(a)); }

void apply_edges_generator::on_input_done() { m_apply_edges.reset(); }

void bind_apply_edges(module& m) {
  define_pipeline_generator<apply_edges_generator, object /* input */,
                            object /* on_assembly_edges */>(m, "apply_edges", "ApplyEdgesGenerator",
                                                            arg("input"), arg("on_assembly_edges"),
                                                            R"DOC(how does this work?)DOC");
}

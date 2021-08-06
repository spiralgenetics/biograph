#include "python/biograph/variants/apply_graph.h"

#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) apply_graph_generator
    : public par_asm_pipeline_wrapper {
 public:
  apply_graph_generator(object input, object on_context);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  object m_on_context;
  std::unique_ptr<apply_graph> m_apply_graph;
};

apply_graph_generator::apply_graph_generator(object input, object on_context)
    : par_asm_pipeline_wrapper(input), m_on_context(on_context) {}

void apply_graph_generator::init() {
  m_apply_graph = make_unique<apply_graph>(
      [this](const graph_context& ctx) {
        gil_scoped_acquire gil;

        m_on_context(ctx);
      },
      make_pipeline_output());
}

void apply_graph_generator::on_input(assembly_ptr a) { m_apply_graph->add(std::move(a)); }

void apply_graph_generator::on_input_done() { m_apply_graph.reset(); }

void bind_apply_graph(module& m) {
  define_pipeline_generator<apply_graph_generator, object /* input */, object /* on_context */>(
      m, "apply_graph", "ApplyGraphGenerator", arg("input"), arg("on_graph_context"),
      R"DOC(how does this work?)DOC");

  class_<graph_context>(m, "GraphContext", "Contains the context of a variant in the variant graph")
      .def_readonly("a", &graph_context::a, "Variant assembly")
      .def_readonly("left_ref", &graph_context::left_ref,
                    "Reference assembly to the left of the variant")
      .def_readonly("refs", &graph_context::refs, "Reference assemblies on the reference branch")
      .def_readonly("right_ref", &graph_context::right_ref,
                    "Reference assemblies to the right of the variant")

      .def("ref_coverage", &graph_context::ref_coverage)
      .def("ref_pair_coverage", &graph_context::ref_coverage)
      .def("ref_scaffold", &graph_context::ref_scaffold)
      .def("edge_coverage",
           static_cast<edge_coverage_t (graph_context::*)(
               const variants::scaffold&, const read_coverage_t&, const read_coverage_t&) const>(
               &graph_context::edge_coverage),
           arg("ref_sequence"), arg("var_cov"), arg("ref_cov"))
      .def(
          "edge_coverage",
          static_cast<edge_coverage_t (graph_context::*)(
              const read_coverage_t&, const read_coverage_t&) const>(&graph_context::edge_coverage),
          arg("var_cov"), arg("ref_cov"));
}

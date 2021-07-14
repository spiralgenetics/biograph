#include "python/biograph/variants/graph_discover.h"

#include "modules/graph_discover/branch.h"
#include "modules/graph_discover/graph_trim_ref.h"
#include "modules/graph_discover/make_ref.h"
#include "modules/graph_discover/push_to_pair.h"
#include "modules/graph_discover/update_rc_seqset_entries.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/par_pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) graph_discover_branch_generator
    : public pipeline_wrapper_with_args<par_asm_pipeline_wrapper, std::shared_ptr<readmap>,
                                        int /* min overlap */, std::string /* tag */,
                                        std::vector<std::string> /* discover_tags */> {
 public:
  using pipeline_wrapper_with_args::pipeline_wrapper_with_args;
  pipeline_step_t make_pipeline_step() override {
    m_options.readmap = std::get<0>(m_args).get();
    m_options.seqset = m_options.readmap->get_seqset().get();
    m_options.min_overlap = std::get<1>(m_args);
    m_options.discover_tags = string_set(std::get<3>(m_args));
    return make_unique<branch_discover>(m_options, std::get<2>(m_args) /* tag */,
                                        this->make_pipeline_output());
  }

 private:
  assemble_options m_options;
};

class __attribute__((visibility("hidden"))) graph_discover_push_to_pair_generator
    : public pipeline_wrapper_with_args<par_asm_pipeline_wrapper, std::shared_ptr<readmap>,
                                        int /* min overlap */, int /* max pair distance */,
                                        std::string /* tag */,
                                        std::vector<std::string> /* discover tags */> {
 public:
  using pipeline_wrapper_with_args::pipeline_wrapper_with_args;
  pipeline_step_t make_pipeline_step() override {
    m_options.readmap = std::get<0>(m_args).get();
    m_options.seqset = m_options.readmap->get_seqset().get();
    m_options.min_overlap = std::get<1>(m_args);
    m_options.max_pair_distance = std::get<2>(m_args);
    m_options.discover_tags = string_set(std::get<4>(m_args));
    return make_unique<push_to_pair_discover>(m_options, std::get<3>(m_args),
                                              this->make_pipeline_output());
  }

 private:
  assemble_options m_options;
};

class __attribute__((visibility("hidden"))) graph_discover_update_rc_seqset_entries_generator
    : public pipeline_wrapper_with_args<par_asm_pipeline_wrapper, std::shared_ptr<seqset>,
                                        bool /* self test */> {
 public:
  using pipeline_wrapper_with_args::pipeline_wrapper_with_args;
  pipeline_step_t make_pipeline_step() override {
    m_options.seqset = std::get<0>(m_args).get();
    auto res = make_unique<update_rc_seqset_entries>(m_options, this->make_pipeline_output());
    if (std::get<1>(m_args)) {
      res->enable_self_test();
    }
    return res;
  }

 private:
  assemble_options m_options;
};
class __attribute__((visibility("hidden"))) graph_discover_trim_ref_generator
    : public pipeline_wrapper_with_args<par_asm_pipeline_wrapper,
                                        const reference_wrapper& /* ref */,
                                        const std::string& /* scaffold_name */, bool /* rev_comp */
                                        > {
 public:
  using pipeline_wrapper_with_args::pipeline_wrapper_with_args;
  pipeline_step_t make_pipeline_step() override {
    m_options.scaffold_name = std::get<1>(m_args);
    m_scaffold = trace_ref::ref_to_scaffold(std::get<0>(m_args).get_reference().get(),
                                            m_options.scaffold_name);
    m_options.scaffold = &m_scaffold;
    if (std::get<2>(m_args)) {
      m_scaffold = m_scaffold.rev_comp();
    }

    return make_unique<graph_trim_ref>(m_options, this->make_pipeline_output());
  }

 private:
  assemble_options m_options;
  ::variants::scaffold m_scaffold;
};

std::vector<assembly_ptr> make_ref_assemblies_wrapper(const reference_wrapper& ref,
                                                      const std::string& scaffold_name,
                                                      aoffset_t start_offset,
                                                      aoffset_t limit_offset,
                                                      aoffset_t max_chunk_size) {
  ::variants::scaffold s = trace_ref::ref_to_scaffold(ref.get_reference().get(), scaffold_name);
  return make_ref_assemblies(s, start_offset, limit_offset, max_chunk_size);
}

void bind_graph_discover(module& m) {
  define_pipeline_generator<graph_discover_branch_generator, object /* input */,
                            const std::shared_ptr<readmap>&, int /* min_overlap */,
                            std::string /* tag */, std::vector<std::string> /* discover tags */>(
      m, "discover_branch", "DiscoverBranchGenerator", arg("input"), arg("readmap"),
      arg("min_overlap") = 70, arg("tag") = "GRAPH_BRANCH",
      arg("discover_tags") = std::vector<std::string>(), R"DOC(how does this work?)DOC");
  define_pipeline_generator<graph_discover_push_to_pair_generator, object /* input */,
                            const std::shared_ptr<readmap>&, int /* min_overlap */,
                            int /* max_pair_distance */, std::string /* tag */,
                            std::vector<std::string> /* discover tags */>(
      m, "discover_push_to_pair", "DiscoverPushToPairGenerator", arg("input"), arg("readmap"),
      arg("min_overlap") = 70, arg("max_pair_distance") = 1000, arg("tag") = "GRAPH_PUSH_TO_PAIR",
      arg("discover_tags") = std::vector<std::string>(), R"DOC(how does this work?)DOC");
  define_pipeline_generator<graph_discover_update_rc_seqset_entries_generator, object /* input */,
                            const std::shared_ptr<seqset>&, bool /* enable self test */>(
      m, "update_rc_seqset_entries", "UpdateRcSeqsetEntriesGenerator", arg("input"), arg("seqset"),
      arg("enable_self_test") = false, R"DOC(how does this work?)DOC");
  define_pipeline_generator<graph_discover_trim_ref_generator, object /* input */,
                            const reference_wrapper& /* ref */,
                            const std::string& /* scaffold_name */, bool /* rev_comp */>(
      m, "graph_trim_ref", "GraphRefTrimGenerator", arg("input"), arg("ref"), arg("scaffold_name"),
      arg("rev_comp") = false, R"DOC(how does this work?)DOC");
  m.def("make_ref_assemblies", make_ref_assemblies_wrapper, arg("ref"), arg("scaffold_name"),
        arg("start_offset"), arg("limit_offset"), arg("max_chunk_size") = 100);
}

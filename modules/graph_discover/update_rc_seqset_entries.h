#pragma once

#include "modules/graph_discover/discover.h"
#include "modules/variants/apply_edges.h"
#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

// Updates rc_seqset_entries in the assemblies
class update_rc_seqset_entries : public apply_edges_step {
 public:
  update_rc_seqset_entries(const assemble_options &options, pipeline_step_t output);
  update_rc_seqset_entries() = delete;
  ~update_rc_seqset_entries() override;

  virtual void on_assembly_edges(optional_aoffset reference_pos,
                                 const std::vector<assembly_ptr> &left_edges,
                                 const std::vector<assembly_ptr> &inserts,
                                 const std::vector<assembly_ptr> &right_edges) override;

  void enable_self_test() { m_self_test = true; }
  bool self_test_succeeded() {
    bool succeeded = m_self_test_succeeded;
    m_self_test_succeeded = true;
    return succeeded;
  }

 private:
  void propagate(const seqset_range_set &incoming, const std::vector<assembly_ptr> &targets);

  assemble_options m_options;
  bool m_self_test = false;
  bool m_self_test_succeeded = true;
};

}  // namespace variants

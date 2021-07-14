// A class to provide vargraph-like coverage on assemblies.
//
// Inputs are:
//   A list of assemblies
//

#include "modules/variants/apply_edges.h"
#include "modules/variants/assemble.h"

namespace variants {

class pair_edge_cov : public apply_edges_step {
 public:
  pair_edge_cov(const assemble_options& opts, pipeline_step_t output);
  ~pair_edge_cov();
  void on_assembly(assembly_ptr a) override;
  void on_advance(aoffset_t new_cur_offset) override;
  void on_assembly_edges(optional_aoffset reference_pos, const std::vector<assembly_ptr>& left_edges,
                         const std::vector<assembly_ptr>& inserts,
                         const std::vector<assembly_ptr>& right_edges) override;

 private:
  static void add_var_edge_coverage(assembly* a);
  static void add_edge_read_ids(const assembly* a, aoffset_t pos, read_id_set* read_ids);

  assemble_options m_options;
};

}  // namespace variants

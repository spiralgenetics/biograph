// A class to provide vargraph-like coverage on assemblies.
//
// Inputs are:
//   A list of assemblies
//

#include "absl/container/btree_map.h"
#include "modules/variants/apply_edges.h"
#include "modules/variants/assemble.h"

namespace variants {

class align_count : public apply_edges_step {
 public:
  align_count(const assemble_options& opts, pipeline_step_t output);
  ~align_count();

  void on_assembly_edges(optional_aoffset reference_pos,
                         const std::vector<assembly_ptr>& left_edges,
                         const std::vector<assembly_ptr>& inserts,
                         const std::vector<assembly_ptr>& right_edges) override;


 private:
  using read_infos_t = absl::btree_map<uint32_t /* read id */, size_t /* count of aligned bases */>;
  struct active_assembly {
    read_id_set all_read_ids;

    align_count_t counts;
  };
  using actives_t = absl::btree_map<assembly*, active_assembly>;

  void start_assembly(assembly* a);
  void end_assembly(assembly* a);

  void add_coverage(assembly* a, active_assembly* act,  bool add /* true to add, false to remove */);

  // Active assemblies
  actives_t m_active;

  // Total of active_assembly::counts for all active assemblies
  absl::btree_map<uint32_t /* read_id */, size_t /* aligned_bases */> m_active_counts;

  assemble_options m_opts;
};

}  // namespace variants

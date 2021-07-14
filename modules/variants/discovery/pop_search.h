#pragma once

#include "modules/variants/discovery/branch.h"
#include "modules/variants/discovery/path.h"

namespace variants {
namespace discovery {

class pop_search_entry : public branch_search_entry {
 public:
  pop_search_entry(const seqset_range& popped, path rc_path, unsigned pair_match_count);

  void check_invariants(const branch* br) const override;
  unsigned cur_overlap() const override { return m_popped.size(); }

  const path& get_path() const override { return m_rc_path; }

 protected:
  search_result search_internal(branch* br) override;
  std::string describe_internal(const branch* br) const override;

 private:
  friend class discovery_test;
  friend class pop_tracer;

  seqset_range m_popped;
  path m_rc_path;
};

}  // namespace discovery
}  // namespace variants

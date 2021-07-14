#pragma once

#include "modules/variants/discovery/branch.h"
#include "modules/variants/discovery/path.h"

namespace variants {
namespace discovery {

class push_search_entry : public branch_search_entry {
 public:
  push_search_entry(path p, unsigned pair_match_count);

  void check_invariants(const branch* br) const override;
  unsigned cur_overlap() const override { return m_path.range().size(); }

  const path& get_path() const override { return m_path; }

 protected:
  search_result search_internal(branch* br) override;
  std::string describe_internal(const branch* br) const override;

 private:
  friend class discovery_test;
  friend class push_tracer;

  void populate_last_read();

  bool check_rejoin();

  // Path starting at br->right_push_view_offset in reference and going to the left.
  path m_path;

  // True if we just did a rejoin.
  bool m_did_rejoin = false;
};

}  // namespace discovery
}  // namespace variants

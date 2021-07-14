#pragma once

#include "modules/variants/discovery/branch.h"
#include "modules/variants/discovery/path.h"

namespace variants {
namespace discovery {

// Rejoins use the push-front view of the branch.
class rejoin_search_entry : public branch_search_entry {
 public:
  static const char k_tracer_name[];
  rejoin_search_entry(unsigned path_overlap, aoffset_t left_offset, aoffset_t left_anchor_len,
                      path p, unsigned pair_match_count);

  void check_invariants(const branch* br) const override;

  const path& get_path() const override { return m_path; }

 protected:
  search_result search_internal(branch* br) override;
  std::string describe_internal(const branch* br) const override;

  void notify_discard(branch* br) override;

 private:
  friend class discovery_test;

  void output(branch* br, bool walk_more_vars);

  aoffset_t const m_left_offset;
  unsigned const m_left_anchor_len;
  path const m_path;
};

}  // namespace discovery
}  // namespace variants

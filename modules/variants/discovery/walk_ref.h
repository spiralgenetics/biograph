#pragma once

#include "modules/variants/discovery/seqset_range_table.h"
#include "modules/variants/discovery/view.h"

namespace variants {
namespace discovery {

class walk_ref_t {
 public:
  walk_ref_t() = delete;
  walk_ref_t(const walk_ref_t&) = delete;
  walk_ref_t(view_t* v) : m_view(v) {}

  void init_pairs_and_push();
  void walk_ref(aoffset_t start, aoffset_t limit);
  void check_invariants() const;

  const assemble_options& opts() const { return m_view->opts(); }

 private:
  struct wr_range_info_t {
    seqset_range r;
    aoffset_t offset;
    boost::optional<dna_base> next_ref_base;
    dna_slice seq;
    bool saved = false;
    int bases_since_read = std::numeric_limits<int>::max() / 2;
  };

  using wr_state_t = std::vector<wr_range_info_t>;

  wr_range_info_t* add_ref_range(aoffset_t pos, const seqset_range& r, dna_slice seq,
                                 boost::optional<dna_base> next_ref_base, int bases_since_read);
  void save_ref_range(wr_range_info_t*);

  std::vector<seqset_range> m_seen_ranges;

  view_t* const m_view;
  wr_state_t m_wr_range_info;
};

}  // namespace discovery
}  // namespace variants

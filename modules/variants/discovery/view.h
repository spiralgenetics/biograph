#pragma once

#include "modules/bio_base/readmap.h"
#include "modules/variants/discovery/seqset_range_table.h"
#include "modules/variants/discovery/state.h"

namespace variants {
namespace discovery {

struct right_partial {
  right_partial() = delete;
  right_partial(dna_slice new_seq, aoffset_t new_outer_right_offset, unsigned new_pair_match_count);

  // Sequence of bases to the right, starting with range.sequence().begin().
  dna_sequence seq;

  // Position in reference that lines up to seq.end().  Anchor starts
  // to the left of right_offset.
  aoffset_t outer_right_offset = std::numeric_limits<aoffset_t>::max();

  // Number of pairs that were found supporting this partially
  // anchored sequence.
  unsigned pair_match_count = std::numeric_limits<unsigned>::max();

  friend std::ostream& operator<<(std::ostream& os, const right_partial& rp) {
    return os << "right_partial@" << rp.outer_right_offset << ":" << rp.seq;
  }
};

// Information associated with a given seqset_range
struct range_info_t {
  // Any offsets where we've seen this seqset_range.  The offset in
  // reference_offsets corresponds to the position of
  // range.sequence().begin().
  std::vector<aoffset_t> reference_offsets;

  // A set of all reference offsets, relative to the direction being
  // processed, where pairing data has indicated that this sequence
  // might be.  The offsets are possible offsets for
  // range.sequence().begin().
  interval_set_t pair_supported_offsets;

  std::vector<right_partial> right_partials;
};

using range_info_table_t = seqset_range_table<range_info_t>;
std::ostream& operator<<(std::ostream& os, const range_info_t& ri);

// A view of the scaffold we're processing.  There are two views: forward ("fwd"), and reverse
// complement ("rev").
class branch;
class path;
class view_t {
 public:
  // Returns pair of fwd view, reverse view.
  static std::pair<std::unique_ptr<view_t>, std::unique_ptr<view_t>> create_view_pair(
      const scaffold& s, state* st);

  bool is_rev_comp() const { return m_is_rev_comp; }
  state* get_state() const { return m_state; }

  void add_push_traces(const path& p, aoffset_t right_offset,
                       boost::optional<dna_base> base_to_skip, branch* on_branch = nullptr);

  void walk_assembly_variants(unsigned path_overlap, aoffset_t left_offset,
                              aoffset_t left_anchor_len, aoffset_t right_offset,
                              aoffset_t right_anchor_len, dna_slice seq, branch* on_branch);

  void add_right_partial(const seqset_range& r, right_partial rp);

  // Get or create a branch off of reference.  The branch extends to the left at right_offset.
  branch* get_branch(dna_base b, aoffset_t right_offset);

  // Same as state::get_offset_info, except relativew to this view.
  offset_info get_offset_info(aoffset_t position, bool fwd) const;

  view_t* reverse_view() const { return m_reverse; }
  aoffset_t reverse_offset(aoffset_t offset) const { return m_scaffold.end_pos() - offset; }
  const scaffold& get_scaffold() const { return m_scaffold; }

  const assemble_options& opts() const { return m_state->opts(); }

  range_info_table_t& range_info() { return m_range_info; }

  // Testing access:
  void check_invariants() const;

  // Returns the number of bases that are shared between the reference
  // to the right of ref_offset and the beginning of slice.
  unsigned shared_ref_bases_to_right(aoffset_t ref_offset, dna_slice slice) const;
  // Returns the number of bases that are shared between the reference
  // to the left of ref_offset and the end of slice.
  unsigned shared_ref_bases_to_left(aoffset_t ref_offset, dna_slice slice) const;

  // Finds all the reads for the given range and adds pair support for
  // all those reads' mates.
  void add_pair_offset_support_for_range(aoffset_t min_offset, aoffset_t max_offset,
                                         const seqset_range& r);

  bool has_ploids_remaining(aoffset_t left_offset, aoffset_t right_offset) const;

  std::vector<branch*> branches() const;

  // Testing access:
  void discard_search_entries();

 private:
  friend class discovery_test;

  view_t() = default;

  void add_ref_range(aoffset_t offset, seqset_range r);
  void add_pair_offset_support(aoffset_t min_offset, aoffset_t max_offset, const readmap::read& rd);

  bool m_is_rev_comp;
  view_t* m_reverse = nullptr;

  scaffold m_scaffold;

  range_info_table_t m_range_info;

  using branch_key = std::pair<dna_base /* base to the left */, aoffset_t /* right offset */>;
  struct branch_hash {
    size_t operator()(const branch_key& k) const {
      std::hash<aoffset_t> offset_hash;
      size_t val = offset_hash(k.second);
      boost::hash_combine(val, int(k.first));
      return val;
    }
  };
  std::unordered_map<branch_key, std::shared_ptr<branch>, branch_hash> m_branches;

  state* m_state = nullptr;
};

}  // namespace discovery
}  // namespace variants

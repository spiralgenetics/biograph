#pragma once

#include <boost/optional.hpp>
#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_set.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/graph_discover/discover.h"
#include "modules/variants/assemble.h"
#include "modules/variants/ref_map.h"
#include "modules/variants/scaffold.h"

namespace variants {

class push_to_pair_discover : public graph_discover {
 public:
  push_to_pair_discover(const assemble_options& options, const std::string& tag,
                        pipeline_step_t output);
  ~push_to_pair_discover() override;

  void flush() override;

 private:
  class entry_anchor_adder;
  struct entry_anchor_search;
  using rc_starts_t = absl::btree_multimap<seqset_range, potential_anchor>;
  struct entry_info {
    rc_starts_t rc_starts;
  };
  struct anchor_info {
    seqset_range truncated;
    seqset_range orig;
    potential_anchor anchor;
  };
  struct push_active_assembly : public active_assembly {
    // Mates supported by this assembly.
    seqset_range_set mates;

    std::vector<anchor_info> rc_anchors;

    std::string to_string() const override;
    ~push_active_assembly() override = default;
  };

  struct potential_entry_anchor : public potential_anchor {
    seqset_range r;

    // Shared bases between the entry searched for and this entry, so
    // we don't have to call .sequence() which is slow.
    aoffset_t shared_bases = 0;

    friend std::ostream& operator<<(std::ostream& os, const potential_entry_anchor& anchor) {
      return os << "r=" << anchor.r.sequence() << " shared=" << anchor.shared_bases << ", "
                << (const potential_anchor&)anchor;
    }
  };

  std::unique_ptr<active_assembly> make_active_assembly() override {
    return make_unique<push_active_assembly>();
  }

  void on_walk(active_assembly*) override;
  void on_readahead(const active_assembly*) override;
  void on_readahead_done(const active_assembly*) override;
  void on_trace(const active_assembly*) override;
  void on_advance_trace(aoffset_t offset) override;

  // Returns all elements from m_potential_rc_anchor_entries that
  // share min_overlap bases with r.
  std::vector<potential_entry_anchor> get_entry_anchors(const seqset_range_set& rs,
                                                        aoffset_t min_overlap,
                                                        aoffset_t approx_ref_offset) const;

  // Returns shared bases
  aoffset_t save_anchors(const active_assembly* act, const seqset_range_set& rs, dna_slice seq,
                         const absl::btree_map<aoffset_t, seqset_range_set>& path_entries,
                         aoffset_t& best_abs_svlen, assembly_ptr& best_abs_svlen_a, aoffset_t& best_shared_bases, assembly_ptr& best_shared_bases_a);

  static seqset_path path_entries_to_rc_path(
      const absl::btree_map<aoffset_t, seqset_range_set>& path_entries, aoffset_t end_pos);

  seqset_range_set trace_one_base(const seqset_range_set& rs, dna_sequence* seq,
                                  int* bases_since_read);

  // Mates we're expecting, and reference count for items in m_mate_expiry
  // that have indicated support for them.
  absl::btree_map<seqset_range, int> m_mates;

  // Offsets and when we remove pair supports m_mates.
  absl::btree_map<aoffset_t, seqset_range_set> m_mate_expiry;

  absl::btree_map<seqset_range /* truncated */, entry_info> m_rc_entry_anchors;

  absl::flat_hash_set<seqset_range, seqset_range_hash> m_seen_ranges;

  // Tag to put on discovered assemblies.
  std::string m_tag;
};

};  // namespace variants

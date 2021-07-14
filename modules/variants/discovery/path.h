#pragma once

#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"

namespace variants {
namespace discovery {

// State of a path being traced.  Calculates overlap on the way.
class path {
 public:
  path() = default;
  path(const path&) = default;
  path(const readmap* rm, dna_slice seq, seqset_range r, unsigned path_overlap,
       unsigned bases_since_read, int anchor_len);
  path& operator=(const path&) = default;

  int bases_since_read() const { return m_bases_since_read; }

  size_t size() const { return m_rc_seq.size(); }

  bool valid() const { return m_r.valid(); }

  bool loop_detected() const { return m_loop_detected; }

  int anchor_len() const { return m_anchor_len; }

  unsigned tot_overlap() const { return m_tot_overlap; }
  unsigned tot_overlap_bases() const { return m_tot_overlap_bases; }

  void push_front_drop(dna_base b);

  // Variant of push_front_drop that allows specifying the pushed
  // range to avoid having to call push_front_drop again.
  void push_front_drop(dna_base b, const seqset_range& pushed);

  // Push the whole sequence on the front.
  void push_front_drop(dna_slice seq);

  int path_overlap() const { return m_path_overlap; }
  int last_overlap() const { return m_last_overlap; }
  int cur_overlap() const { return int(m_r.size()) - bases_since_read(); }
  const seqset_range& range() const { return m_r; }
  dna_slice seq() const { return dna_slice(m_rc_seq).rev_comp(); }
  const boost::optional<uint32_t> longest_read_id() const { return m_longest_read_id; }

  friend std::ostream& operator<<(std::ostream& os, const path& p) {
    p.display_fwd(os);
    return os;
  }

  void display_fwd(std::ostream& os) const;
  void display_rev(std::ostream& os) const;

  // For testing:
  void set_path_overlap(int new_path_overlap) { m_path_overlap = new_path_overlap; }

  void check_invariants() const;

 private:
  void populate_longest_read();

  const readmap* m_rm = nullptr;
  dna_sequence m_rc_seq;
  seqset_range m_r;
  int m_path_overlap = 0;
  int m_last_overlap = 0;
  int m_bases_since_read = 0;
  int m_anchor_len = 0;

  // sum(cur_overlap) over tot_overlap_bases
  unsigned m_tot_overlap = 0;
  unsigned m_tot_overlap_bases = 0;

  boost::optional<uint32_t> m_longest_read_id;

  unsigned m_bases_until_loop_check = 1;
  unsigned m_loop_check_interval = 1;
  uint64_t m_loop_check_seqset_id = std::numeric_limits<uint64_t>::max();
  bool m_loop_detected = false;
};

}  // namespace discovery
}  // namespace variants

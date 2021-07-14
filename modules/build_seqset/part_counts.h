#pragma once

#include "modules/build_seqset/repo_seq.h"

namespace build_seqset {

// Counts distribution statistics for beginnings of sequences.
class part_counts {
 public:
  // Number of bases to count individually.
  part_counts(unsigned n_bases);
  part_counts() = delete;
  part_counts(const part_counts&) = delete;

  std::pair<size_t, size_t> seq_to_index_range(dna_slice seq);

  size_t get_index(const seq_repository::entry_data& e) const {
    auto it = dna_const_iterator(e.raw_inline_bases(), 0, 0 /* rev comp */);
    size_t idx = 0;
    for (unsigned i = 0; i != m_bases; ++i) {
      idx <<= 2;
      idx |= int(dna_base(*it));
      ++it;
    }
    CHECK_LT(idx, m_counts.size()) << e;
    return idx;
  }

  void add(const seq_repository::entry_data& e) {
    size_t idx = get_index(e);
    ++m_counts[idx];
  }

  const tracked_vector<size_t>& counts() const { return m_counts; }
  size_t bases() const { return m_bases; }

  std::string display_histo() const;

 private:
  unsigned const m_bases;

  tracked_vector<size_t> m_counts;
};

}  // namespace build_seqset

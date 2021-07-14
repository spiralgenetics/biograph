#pragma once

#include <queue>
#include <vector>
#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_base/seqset_mergemap.h"
#include "modules/io/progress.h"

// make_mergemap counts shared sequences between multiple input flat
// seqsets.
//
// Let's call the input part seqsets P1, P2, P3, ....
// We want to merge them into "M", the resultant merged seqset (which hasn't
// been constructed yet).
//
// make_mergemap constructs bitcounts B1, B2, B3, ... for each P1,
// P2, P3 ...
//
// Each of these bitcounts will be of the size M.
//
// A bit x will be set in B1 if the sequence with index x in M (or
// something beginning with it; see below) is present in P1.  In this
// case, the sequence with index x in M will be the sequence with
// index B1.count(x) in P1.
//
// Prefixes:
//
// Seqsets guarantee that if sequence X exists as an entry in the
// seqset, there is no other sequence that starts with X that also
// exists as an entry in the seqset.
//
// This means that if X is a prefix of Y, and X exists in P1 and Y
// exists in P2, that only Y can exist in M.
//
// In this case, the bit corresponding to X will be set in B1, and the
// bit corresponding to Y will be set in B2.
//
// Implementation:
//
// In order to parallelize, make_mergemap must split up
// the input seqsets.  It does this by selecting the largest input part in
// terms of number of sequences, and splitting based on sequence index.
//
// It then does a binary search on each of the smaller input parts in
// order to find the same sequence range from the input.

class make_mergemap {
 public:
  make_mergemap(const std::vector<const seqset_flat*>& flats);

  void build(progress_handler_t progress = null_progress_handler);

  size_t total_merged_entries() const { return m_total_merged_entries; }
  void fill_mergemap(unsigned input_id, seqset_mergemap_builder* mergemap,
                     progress_handler_t progress = null_progress_handler);

  static int64_t queue_entry_compare(const dna_slice& lhs, unsigned lhs_size,
                                     const dna_slice& rhs, unsigned rhs_size);
  // Entry in priority queue for sorted merging.
  struct queue_entry {
    dna_slice cur_slice;
    unsigned flat;  // index into m_flat
    size_t entry_id;

    bool operator<(const queue_entry& rhs) const {
      // The top of the prioriry queue should be the smallest element,
      // not the largest, so reverse the comparison:
      return cur_slice > rhs.cur_slice;
    }
  };

  void count_range(size_t start, size_t limit);

  enum class empty_meaning { start, limit };
  friend std::ostream& operator<<(std::ostream& os, empty_meaning empty_means);

  dna_slice seq_for_pos(size_t pos, empty_meaning empty_means) const;
  void positions_for_seq(dna_slice slice, std::vector<size_t>& indees,
                         empty_meaning empty_means) const;

  void add_to_queue(std::priority_queue<queue_entry>& queue, unsigned flat,
                    size_t entry_id, size_t limit_idx) const;

 public:
  // Result from an input part
  struct chunk_result {
    // Entry count of merged entries
    size_t merged_entries = 0;

    // Results corresponding to entries in m_flats.
    std::vector<std::vector<bool>> bits;
  };

  std::vector<const seqset_flat*> m_flats;

 private:
  // Input flat with the largest size.
  unsigned m_biggest_flat;

  std::mutex m_mu;
  std::map<dna_slice, chunk_result> m_chunk_results;
  size_t m_total_merged_entries = 0;
};

inline std::ostream& operator<<(std::ostream& os,
                                make_mergemap::empty_meaning empty_means) {
  switch (empty_means) {
    case make_mergemap::empty_meaning::start:
      return os << "start";
    case make_mergemap::empty_meaning::limit:
      return os << "limit";
  }
  LOG(FATAL) << int(empty_means);
  return os;
}

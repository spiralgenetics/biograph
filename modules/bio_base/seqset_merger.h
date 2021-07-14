#pragma once

#include <vector>

#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_base/seqset_mergemap.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file.h"

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

class seqset_merger {
 public:
  seqset_merger(const std::vector<const seqset_flat*>& flats,
                const std::vector<const seqset_mergemap*>& mergemaps);

  void build(const spiral_file_create_state& state,
             progress_handler_t progress = null_progress_handler);

  class iterator
      : public boost::iterator_facade<iterator, dna_slice const,
                                      std::random_access_iterator_tag,
                                      dna_slice, ptrdiff_t> {
   public:
    iterator() : m_merger(nullptr), m_merge_idx(0) {}
    iterator(const seqset_merger* merger, size_t merge_idx)
        : m_merger(merger),
          m_merge_idx(merge_idx),
          m_part_idx(merger->m_flats.size()) {
      seek_to_merge_idx();
    }
    bool equal(const iterator& rhs) const {
      return m_merge_idx == rhs.m_merge_idx;
    }
    void increment();
    void decrement();
    void advance(ptrdiff_t offset) {
      m_merge_idx += offset;
      seek_to_merge_idx();
    }
    ptrdiff_t distance_to(const iterator& rhs) const {
      return rhs.m_merge_idx - m_merge_idx;
    }
    dna_slice dereference() const;

    size_t merge_idx() const { return m_merge_idx; }

   private:
    void seek_to_merge_idx();
    const seqset_merger* m_merger;
    size_t m_merge_idx;
    std::vector<size_t> m_part_idx;
  };

  iterator begin() const { return iterator(this, 0); }
  iterator end() const { return iterator(this, m_num_seqs); }

 private:
  void merge_range(size_t start, size_t limit);
  iterator get_base_iterator(dna_base base, iterator it);

  // Inputs:
  std::vector<const seqset_flat*> m_flats;
  std::vector<const seqset_mergemap*> m_mergemaps;
  unsigned m_num_inputs = 0;
  size_t m_num_seqs = 0;

  // Output:
  std::unique_ptr<seqset> m_seqset;
};

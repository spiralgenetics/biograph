#pragma once

#include "modules/bio_base/readmap.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/bio_base/seqset.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file.h"
#include "modules/variants/ref_map.h"

namespace variants {

// Samples pair data statistics from a seqset.  Samples a bunch of
// reads where both ends uniquely match reference, and calculates an
// average pair distance.
//
// There are two different read strategies for dealing with pairs:
//
// "paired end reads", where the reads of the pair both point inwards
// "mate pair reads", where the reads of the pair both point outwards
//
// For more information see:
// https://www.illumina.com/documents/products/datasheets/datasheet_genomic_sequence.pdf
//
// For "paired end reads", the pair distance is the distance from the
// beginning of a forward read to the beginning of its mate.  This is
// the end-to-end distance encompassing a read, its mate, and all the
// bases that span them.
//
// For "mate pair reads", the pair distance is the negative distance
// of the bases that span in between the reads.
//
// This means that in either case, if a read points forward in the
// reference, the offset in the reference of the beginning of its mate
// is the offset of the read plus the pair distance.
//
// If no pairs are found, the average pair distance is 0.
class pair_stats {
 public:
  pair_stats(const seqset* the_seqset, const readmap* the_readmap,
             const reference* ref, const ref_map* rmap)
      : m_seqset(the_seqset),
        m_readmap(the_readmap),
        m_ref(ref),
        m_ref_map(rmap) {}

  void calc_stats(progress_handler_t progress = null_progress_handler);

  bool found_pairs() const { return m_median_pair_offset.is_initialized(); }

  int64_t median_pair_offset() const {
    CHECK(m_median_pair_offset);
    return *m_median_pair_offset;
  }

 private:
  static constexpr size_t k_num_samples = 1000;
  static constexpr size_t k_max_attempts = k_num_samples * 20;

  void get_ref_loc(uint64_t seqset_id, const ref_map::entry& rmap_entry,
                   seq_position& loc, bool& loc_rc);

  const seqset* m_seqset = nullptr;
  const readmap* m_readmap = nullptr;
  const reference* m_ref = nullptr;
  const ref_map* m_ref_map = nullptr;

  boost::optional<int64_t> m_median_pair_offset;
};

}  // namespace

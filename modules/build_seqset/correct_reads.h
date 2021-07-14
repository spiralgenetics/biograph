#pragma once

#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_mapred/correct_reads_mapper.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/bio_mapred/read_correction.h"
#include "modules/build_seqset/part_repo.h"
#include "modules/io/track_mem.h"
#include "modules/build_seqset/part_counts.h"
#include "modules/bio_base/fast_read_correct.h"

namespace build_seqset {

class correct_reads {
 public:
  correct_reads(part_repo& entries, kmer_set& ks,
                const read_correction_params& params);
  ~correct_reads();

  void add_initial_repo(progress_handler_t progress = null_progress_handler);
  // Returns true if the read was successfully corrected.
  bool correct(const unaligned_read& r, corrected_read& cr);

 private:
  bool check_initial_repo_match(const dna_slice& seq, size_t offset,
                                size_t reference_offset);
  void find_reads_and_reference(const dna_slice& seq, const std::vector<frc_kmer>& kmers,
                                size_t& ref_pos, bool& ref_is_rc);

  bool kmer_starts_read(const frc_kmer& k) const;

  part_repo& m_entries;
  kmer_set& m_ks;
  dna_slice m_initial_repo;
  read_correction_params m_params;

  static constexpr uint32_t k_ref_offset_ambiguous = std::numeric_limits<uint32_t>::max();
  static constexpr uint32_t k_ref_offset_not_present = std::numeric_limits<uint32_t>::max() - 1;
  mutable_packed_vector<uint32_t, sizeof(uint32_t) * 8> m_ref_offsets;
  unsigned m_kmer_size;

  std::atomic<uint64_t> m_reference_match{0};
  std::atomic<uint64_t> m_non_reference_match{0};
  std::atomic<uint64_t> m_false_reference_match{0};
  std::atomic<uint64_t> m_corrected_bases{0};
  std::atomic<uint64_t> m_reads_modified{0};
  std::atomic<uint64_t> m_dropped_bases{0};
  std::atomic<uint64_t> m_reads_truncated{0};

  read_correction_stats m_stats;
};

}  // namespace build_seqset

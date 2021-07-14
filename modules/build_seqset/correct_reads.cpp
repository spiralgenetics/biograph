#include "modules/build_seqset/correct_reads.h"
#include <boost/range/adaptor/reversed.hpp>
#include "modules/bio_base/fast_read_correct.h"
#include "modules/build_seqset/repo_seq.h"
#include "modules/io/parallel.h"

namespace build_seqset {

namespace {

kmer_t kmer_shift_left(kmer_t orig, unsigned kmer_size, dna_base b) {
  kmer_t result = orig;
  result <<= 2;
  result |= int(b);
  result &= ~(std::numeric_limits<kmer_t>::max() << (kmer_size * 2));
  return result;
}

}  // namespace

constexpr uint32_t correct_reads::k_ref_offset_not_present;
constexpr uint32_t correct_reads::k_ref_offset_ambiguous;

correct_reads::correct_reads(part_repo& entries, kmer_set& ks, const read_correction_params& params)
    : m_entries(entries),
      m_ks(ks),
      m_initial_repo(entries.repo_slice()),
      m_params(params),
      m_ref_offsets(ks.size(), "correct_reads:ref_offsets"),
      m_kmer_size(m_ks.kmer_size()) {
  if (m_params.skip_snps) {
    throw(io_exception("skip_snps not supported"));
  }
  if (m_params.trim != 0) {
    throw(io_exception("Trim not supported."));
  }
  if (!m_params.frc_max_corrections && !m_params.exact) {
    throw(io_exception("Slow read correction not supported."));
  }
  CHECK_EQ(m_kmer_size, m_ks.kmer_size());
  for (size_t i = 0; i != m_ks.size(); ++i) {
    m_ref_offsets[i] = k_ref_offset_not_present;
  }
}

correct_reads::~correct_reads() {
  size_t ref = m_reference_match.load();
  size_t nonref = m_non_reference_match.load();
  size_t tot = ref + nonref;
  size_t false_ref = m_false_reference_match.load();
  size_t modified = m_reads_modified.load();
  size_t modified_bases = m_corrected_bases.load();
  size_t truncated = m_reads_truncated.load();
  size_t truncated_bases = m_dropped_bases.load();

  SPLOG(
      "%ld corrected reads processed; %ld (%.2f%%) matched reference, %ld "
      "(%.2f%%) did not, including %ld (%.2f%%) which included a kmer matching "
      "reference but did not entirely match.",
      tot, ref, ref * 100. / tot, nonref, nonref * 100. / tot, false_ref, false_ref * 100. / tot);

  SPLOG(
      "%ld bases were corrected in %ld (%.2f%%) reads, averaging %.2f bases "
      "per corrected read.",
      modified_bases, modified, modified * 1. / tot, modified_bases * 1. / modified);

  SPLOG(
      "%ld bases dropped from the end of %ld (%.2f%%) reads, averaging %.2f "
      "bases per truncated read.",
      truncated_bases, truncated, truncated * 100. / tot, truncated_bases * 1. / truncated);
};

void correct_reads::add_initial_repo(progress_handler_t progress) {
  parallel_for(  //
      0, m_initial_repo.size(),
      [this](size_t start, size_t limit) {
        if (start > (m_kmer_size - 1)) {
          start -= (m_kmer_size - 1);
        } else {
          start = 0;
        }
        if (limit - start < m_kmer_size) {
          return;
        }

        unsigned initial_kmer_left = m_kmer_size;
        auto rit_end = m_initial_repo.begin() + limit;
        kmer_t k = 0;
        auto rit = m_initial_repo.begin() + start;
        while (rit != rit_end) {
          k = kmer_shift_left(k, m_kmer_size, *rit);
          ++rit;

          if (initial_kmer_left) {
            initial_kmer_left--;
          }
          if (!initial_kmer_left) {
            bool flipped;
            canonicalize(k, m_kmer_size, flipped);
            if (!flipped) {
              auto it = m_ks.find(k);
              if (it != m_ks.end()) {
                uint32_t kmer_id = it - m_ks.begin();
                uint32_t offset = rit - m_initial_repo.begin() - m_kmer_size;

                for (;;) {
                  uint32_t orig_offset = m_ref_offsets[kmer_id];
                  uint32_t new_offset = offset;

                  if (orig_offset == new_offset || orig_offset == k_ref_offset_ambiguous) {
                    // Original ok.
                    break;
                  }
                  if (orig_offset != k_ref_offset_not_present) {
                    // Different offset there; mark as ambiguous
                    new_offset = k_ref_offset_ambiguous;
                  }

                  if (m_ref_offsets[kmer_id].compare_and_swap(orig_offset, new_offset)) {
                    // Success!
                    break;
                  }

                  // Otherwise, try again.
                }
              }
            }
          }
        }
      },
      progress);
  size_t kmer_matches_ref = 0;
  size_t ambig_ref = 0;
  for (size_t i = 0; i != m_ref_offsets.size(); ++i) {
    switch (m_ref_offsets[i]) {
      case k_ref_offset_not_present:
        break;
      case k_ref_offset_ambiguous:
        ++ambig_ref;
        break;
      default:
        ++kmer_matches_ref;
        break;
    }
  }
  SPLOG(
      "%ld initial bases present for seqset build.  %ld/%ld kmers (%.2f%%) "
      "matched; %ld (%.2f%%) more match to more than one reference location.",
      m_initial_repo.size(), kmer_matches_ref, m_ref_offsets.size(),
      kmer_matches_ref * 100. / m_ref_offsets.size(), ambig_ref,
      ambig_ref * 100. / m_ref_offsets.size());
}

bool correct_reads::correct(const unaligned_read& r, corrected_read& cr) {
  if (r.sequence.size() < m_kmer_size) {
    return false;
  }

  frc_params params;
  params.max_corrections = m_params.frc_max_corrections;
  params.min_good_run = m_params.frc_min_good_run;
  params.kmer_size = m_kmer_size;
  const kmer_set* ks = &m_ks;
  params.kmer_lookup_f = [ks](kmer_t kmer, frc_kmer* ki) -> bool {
    kmer_t canon = canonicalize(kmer, ks->kmer_size(), ki->flipped);
    auto index = ks->find_table_index(canon);
    if (index == kmer_set::k_not_present) {
      return false;
    }
    ki->index = index;
    return true;
  };
  frc_output res = fast_read_correct(r.sequence, params);
  unsigned needed_good_bases = m_params.trim_after_portion * r.sequence.size();
  if (res.corrected.size() < needed_good_bases) {
    m_stats.failed_correction_count++;
    return false;
  }

  if (res.corrected.size() < r.sequence.size()) {
    m_dropped_bases += r.sequence.size() - res.corrected.size();
    m_reads_truncated++;
  }
  if (res.corrections) {
    m_reads_modified++;
    m_corrected_bases += res.corrections;
  }

  const auto& seq = res.corrected;

  bool ref_is_rc;
  size_t ref_pos = k_ref_offset_not_present;
  find_reads_and_reference(seq, res.kmers, ref_pos, ref_is_rc);

  unsigned next_fwd_read = 0;
  for (const frc_kmer& k : res.kmers) {
    if (next_fwd_read > 0 && kmer_starts_read(k)) {
      break;
    }
    ++next_fwd_read;
  }

  unsigned next_rev_read = 0;
  for (const frc_kmer& fwd_k : boost::adaptors::reverse(res.kmers)) {
    frc_kmer k = fwd_k.as_flipped();
    if (next_rev_read > 0 && kmer_starts_read(k)) {
      break;
    }
    ++next_rev_read;
  }

  CHECK_GT(next_fwd_read, 0);
  CHECK_GT(next_rev_read, 0);

  if (ref_pos < k_ref_offset_not_present) {
    m_reference_match++;
    dna_slice c = res.corrected;
    if (ref_is_rc) {
      c = c.rev_comp();
      std::swap(next_fwd_read, next_rev_read);
    }
    m_entries.write_using_repo(c, next_fwd_read, next_rev_read, ref_pos);
  } else {
    m_non_reference_match++;
    m_entries.write(res.corrected, next_fwd_read, next_rev_read);
  }
  cr.corrected = res.corrected;
  m_stats.corrected_read_count++;
  m_stats.corrected_read_bases += res.corrected.size();
  return true;
}

void correct_reads::find_reads_and_reference(const dna_slice& seq,
                                             const std::vector<frc_kmer>& kmers, size_t& ref_pos,
                                             bool& ref_is_rc) {
  ref_pos = k_ref_offset_not_present;

  bool tried_flipped = false;
  bool tried_unflipped = false;

  unsigned offset = 0;
  for (auto it = kmers.begin(); it != kmers.end(); ++it, ++offset) {
    const frc_kmer& k = *it;

    if (tried_flipped && tried_unflipped) {
      return;
    }

    if (k.flipped) {
      if (tried_flipped) {
        continue;
      }
    } else {
      if (tried_unflipped) {
        continue;
      }
    }

    uint32_t ref_offset_candidate = m_ref_offsets[k.index];

    if (ref_offset_candidate == k_ref_offset_not_present) {
      if (k.flipped) {
        tried_flipped = true;
      } else {
        tried_unflipped = true;
      }

      continue;
    }

    if (ref_offset_candidate == k_ref_offset_ambiguous) {
      continue;
    }

    if (k.flipped) {
      size_t flipped_offset = seq.size() - offset - m_kmer_size;
      if (check_initial_repo_match(seq.rev_comp(), flipped_offset, ref_offset_candidate)) {
        ref_pos = ref_offset_candidate - flipped_offset;
        ref_is_rc = true;
      } else {
        ++m_false_reference_match;
      }
    } else {
      if (check_initial_repo_match(seq, offset, ref_offset_candidate)) {
        ref_pos = ref_offset_candidate - offset;
        ref_is_rc = false;
      } else {
        ++m_false_reference_match;
      }
    }
    return;
  }
}

bool correct_reads::check_initial_repo_match(const dna_slice& seq, size_t offset,
                                             size_t reference_offset) {
  DCHECK_LE(reference_offset + m_kmer_size, m_initial_repo.size());
  if (reference_offset < offset) {
    return false;
  }
  DCHECK_EQ(seq.subseq(offset, m_kmer_size), m_initial_repo.subseq(reference_offset, m_kmer_size));
  if (reference_offset - offset + seq.size() > m_initial_repo.size()) {
    return false;
  }
  return seq == m_initial_repo.subseq(reference_offset - offset, seq.size());
}

bool correct_reads::kmer_starts_read(const frc_kmer& k) const {
  unsigned check_flag = k.flipped ? kmer_set::k_rev_starts_read : kmer_set::k_fwd_starts_read;
  return m_ks.get_flags(k.index) & check_flag;
}

}  // namespace build_seqset

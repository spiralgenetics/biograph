#include "modules/bio_base/seqset_merger.h"
#include "base/base.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_mergemap.h"
#include "modules/io/file_io.h"
#include "modules/io/parallel.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/test/coverage.h"

#include <sys/mman.h>

#include <boost/filesystem.hpp>

DECLARE_TEST_COVERAGE(seqset_merger);

seqset_merger::seqset_merger(
    const std::vector<const seqset_flat*>& flats,
    const std::vector<const seqset_mergemap*>& seqset_mergemaps)
    : m_flats(flats),
      m_mergemaps(seqset_mergemaps),
      m_num_inputs(m_flats.size()) {
  CHECK_EQ(m_flats.size(), m_mergemaps.size());
  CHECK_GT(m_num_inputs, 0);
  m_num_seqs = m_mergemaps[0]->get_bitcount().size();
  for (unsigned i = 0; i < m_num_inputs; ++i) {
    CHECK_EQ(m_flats[i]->get_seqset()->uuid(),
             m_mergemaps[i]->metadata().orig_seqset_uuid);
    CHECK_EQ(m_num_seqs, m_mergemaps[i]->get_bitcount().size());
    CHECK_EQ(m_mergemaps[i]->get_bitcount().total_bits(), m_flats[i]->size());
  }
}

namespace {

unsigned shared_prefix_length(const dna_slice& seq1, const dna_slice& seq2) {
  unsigned shared_length = std::min(seq1.size(), seq2.size());

  for (unsigned i = 0; i < shared_length; ++i) {
    if (seq1[i] != seq2[i]) {
      return i;
    }
  }
  return shared_length;
}

}  // namespace

void seqset_merger::build(const spiral_file_create_state& state,
                          progress_handler_t progress) {
  CHECK(!m_seqset);

  for (unsigned i = 0; i < m_num_inputs; ++i) {
    CHECK_EQ(state.uuid(), m_mergemaps[i]->metadata().merged_seqset_uuid);
  }

  unsigned max_read_len = 0;
  for (const auto& f : m_flats) {
    max_read_len = std::max<unsigned>(max_read_len, f->get_seqset()->max_read_len());
  }

  SPLOG("Creating merged seqset with %ld entries, %d maximum entry length", m_num_seqs, max_read_len);
  m_seqset.reset(new seqset(state, m_num_seqs, max_read_len));

  parallel_for(
      0, m_num_seqs,
      [this](size_t start, size_t limit) { merge_range(start, limit); },
      progress);
  SPLOG("Finalizing merged seqset");
  m_seqset->finalize(progress);
  SPLOG("Done creating merged seqest");
}

seqset_merger::iterator seqset_merger::get_base_iterator(dna_base base,
                                                         iterator it) {
  if (it == end()) {
    if (base == dna_base(3)) {
      NOTE_TEST_COVERAGE(seqset_merger);
      return end();
    }

    NOTE_TEST_COVERAGE(seqset_merger);
    dna_sequence search_for;
    search_for.push_back(dna_base(int(base) + 1));
    return std::lower_bound(begin(), end(),
                            dna_slice(search_for.begin(), search_for.end()));
  }

  dna_sequence search_for;
  search_for.push_back(base);
  search_for += dna_sequence(it->begin(), it->end());
  dna_slice search_for_slice(search_for.begin(), search_for.end());
  iterator result = std::lower_bound(begin(), end(), search_for_slice);
  while (result != begin()) {
    iterator prev = result;
    --prev;
    unsigned compare_len = std::min(prev->size(), search_for.size());
    if (subseq_compare(prev->begin(), search_for.begin(), compare_len,
                       compare_len) == dna_compare_result::EQUAL) {
      NOTE_TEST_COVERAGE(seqset_merger);
      result = prev;
    } else {
      NOTE_TEST_COVERAGE(seqset_merger);
      return result;
    }
  }
  return result;
}

void seqset_merger::merge_range(size_t start_offset, size_t limit_offset) {
  iterator start = begin() + start_offset;
  iterator limit = begin() + limit_offset;

  iterator end_iterator = end();

  dna_slice prev_seq;
  if (start == begin()) {
    NOTE_TEST_COVERAGE(seqset_merger);
  } else {
    NOTE_TEST_COVERAGE(seqset_merger);
    prev_seq = *(start - 1);
  }
  dna_base_array<iterator> base_iterator;
  dna_base_array<iterator> base_limit;
  dna_base_array<dna_slice> base_slice;
  for (dna_base base : dna_bases()) {
    base_iterator[base] = get_base_iterator(base, start);
    base_limit[base] = get_base_iterator(base, limit);
    if (base_iterator[base] == end_iterator) {
      NOTE_TEST_COVERAGE(seqset_merger);
    } else {
      NOTE_TEST_COVERAGE(seqset_merger);
      base_slice[base] = *base_iterator[base];
    }
  }

  for (iterator cur = start; cur != limit; ++cur) {
    dna_slice cur_seq = *cur;
    DCHECK(cur_seq > prev_seq) << "\nCur:  " << cur_seq.as_string()
                               << "\nPrev: " << prev_seq.as_string() << "\n";
    DCHECK_GT(cur_seq.size(), 0);

    m_seqset->set_entry_size(cur.merge_idx(), cur_seq.size());

    for (dna_base base : dna_bases()) {
      // Compare the main DNA sequence with an overlap candidate.  To
      // qualify, the candidate suffix (i.e. everything but the first
      // base) must be a prefix of the main sequence.  If there's an
      // overlap, advance the candidate sequence buffer and return
      // true.
      if (base_iterator[base] == base_limit[base]) {
        NOTE_TEST_COVERAGE(seqset_merger);
        continue;
      }
      unsigned overlap_size =
          std::min(base_slice[base].size() - 1, cur_seq.size());
      // The potential overlap is the smaller of the candidate suffix
      // or main sequence
      DCHECK_EQ(base_slice[base][0], base);

      dna_compare_result compare =
          subseq_compare(base_slice[base].begin() + 1, cur_seq.begin(),
                         overlap_size, overlap_size);
      if (compare == dna_compare_result::EQUAL) {
        m_seqset->set_bit(cur.merge_idx(), base, true);
        base_iterator[base]++;
        if (base_iterator[base] == base_limit[base]) {
          NOTE_TEST_COVERAGE(seqset_merger);
        } else {
          NOTE_TEST_COVERAGE(seqset_merger);
          base_slice[base] = *base_iterator[base];
        }
      } else {
        CHECK_EQ(compare, dna_compare_result::SECOND_IS_LESS)
            << "Out-of-order prevs for overlap " << overlap_size
            << ": \nCur:  " << cur_seq.as_string() << "\n"
            << base << ":   " << base_slice[base].as_string();
        NOTE_TEST_COVERAGE(seqset_merger);
      }
    }

    m_seqset->set_shared(cur.merge_idx(),
                         shared_prefix_length(cur_seq, prev_seq));

    prev_seq = cur_seq;
  }

  for (dna_base base : dna_bases()) {
    CHECK(base_iterator[base] == base_limit[base])
        << base << " " << base_iterator[base].merge_idx() << " "
        << base_limit[base].merge_idx() << "\nM:  " << prev_seq.as_string()
        << "\nL:  " << limit->as_string()
        << "\nb: " << base_iterator[base]->as_string()
        << "\nl: " << base_limit[base]->as_string() << "\n";
  }
}

void seqset_merger::iterator::seek_to_merge_idx() {
  for (unsigned flat_id = 0; flat_id < m_merger->m_flats.size(); flat_id++) {
    m_part_idx[flat_id] =
        m_merger->m_mergemaps[flat_id]->get_bitcount().count(m_merge_idx);
  }
}

void seqset_merger::iterator::increment() {
  for (unsigned flat_id = 0; flat_id < m_merger->m_flats.size(); flat_id++) {
    if (m_merger->m_mergemaps[flat_id]->get_bitcount().get(m_merge_idx)) {
      m_part_idx[flat_id]++;
    }
  }
  m_merge_idx++;
}

void seqset_merger::iterator::decrement() {
  DCHECK_GT(m_merge_idx, 0);
  m_merge_idx--;
  for (unsigned flat_id = 0; flat_id < m_merger->m_flats.size(); flat_id++) {
    if (m_merger->m_mergemaps[flat_id]->get_bitcount().get(m_merge_idx)) {
      m_part_idx[flat_id]--;
    }
  }
}

dna_slice seqset_merger::iterator::dereference() const {
  dna_slice seq;
  for (unsigned flat_id = 0; flat_id < m_merger->m_flats.size(); ++flat_id) {
    if (m_merger->m_mergemaps[flat_id]->get_bitcount().get(m_merge_idx)) {
      dna_slice part_seq = m_merger->m_flats[flat_id]->get(m_part_idx[flat_id]);
      if (part_seq.size() > seq.size()) {
        NOTE_TEST_COVERAGE(seqset_merger);
        seq = part_seq;
      } else {
        NOTE_TEST_COVERAGE(seqset_merger);
      }
    }
  }
  DCHECK_GT(seq.size(), 0);
  return seq;
}

#include "modules/bio_base/make_mergemap.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/io/parallel.h"
#include "modules/test/coverage.h"

#include <bitset>

DECLARE_TEST_COVERAGE(make_mergemap);

namespace {

bool is_equal_or_prefix(dna_slice prefix, dna_slice longer) {
  if (prefix.size() > longer.size()) {
    return false;
  }

  return prefix == longer.subseq(0, prefix.size());
}

}  // namespace

make_mergemap::make_mergemap(const std::vector<const seqset_flat*>& flats)
    : m_flats(flats) {}

void make_mergemap::build(progress_handler_t progress) {
  SPLOG("Creating new seqset flat merge counter for %ld parts", m_flats.size());
  NOTE_TEST_COVERAGE_IF(make_mergemap, m_flats.size() == 1);
  NOTE_TEST_COVERAGE_IF(make_mergemap, m_flats.size() == 2);
  NOTE_TEST_COVERAGE_IF(make_mergemap, m_flats.size() > 2);

  CHECK(!m_flats.empty());
  m_biggest_flat = 0;

  for (unsigned i = 0; i < m_flats.size(); i++) {
    if (m_flats[i]->size() > m_flats[m_biggest_flat]->size()) {
      m_biggest_flat = i;
    }
  }

  parallel_for(
      0, m_flats[m_biggest_flat]->size(),
      [this](size_t start, size_t limit) { count_range(start, limit); },
      progress);
  SPLOG("Done counting parts");
}

dna_slice make_mergemap::seq_for_pos(size_t pos,
                                     empty_meaning empty_means) const {
  if (pos == 0) {
    NOTE_TEST_COVERAGE(make_mergemap);
    CHECK_EQ(empty_meaning::start, empty_means);
    return dna_slice();
  }

  const seqset_flat* flat = m_flats[m_biggest_flat];

  if (pos == flat->size()) {
    NOTE_TEST_COVERAGE(make_mergemap);
    CHECK_EQ(empty_meaning::limit, empty_means);
    return dna_slice();
  }

  NOTE_TEST_COVERAGE(make_mergemap);
  CHECK_LT(pos, flat->size());
  return flat->get(pos);
}

void make_mergemap::positions_for_seq(dna_slice target_seq,
                                      std::vector<size_t>& indexes,
                                      empty_meaning empty_means) const {
  if (target_seq.size() == 0) {
    for (unsigned input_id = 0; input_id < m_flats.size(); ++input_id) {
      const seqset_flat* flat = m_flats[input_id];
      switch (empty_means) {
        case empty_meaning::start:
          NOTE_TEST_COVERAGE(make_mergemap);
          indexes[input_id] = 0;
          break;
        case empty_meaning::limit:
          NOTE_TEST_COVERAGE(make_mergemap);
          indexes[input_id] = flat->size();
          break;
      }
    }
    return;
  }

  // Find the shortest prefix we can find to start the block.  If we find a
  // shorter one, repeat the process.
  int found_shorter_count = 0;
  bool found_shorter = true;
  while (found_shorter) {
    found_shorter = false;

    for (unsigned input_id = 0; input_id < m_flats.size(); ++input_id) {
      const seqset_flat* flat = m_flats[input_id];
      seqset_flat::iterator it =
          std::upper_bound(flat->begin(), flat->end(), target_seq);
      size_t pos = it - flat->begin();
      if (pos != 0) {
        // Check to see if this entry exists in this seqset; if so,
        // start before it.
        size_t maybe_pos = pos - 1;
        dna_slice maybe_prefix = flat->get(maybe_pos);
        if (is_equal_or_prefix(maybe_prefix, target_seq)) {
          if (maybe_prefix.size() < target_seq.size()) {
            // Found a prefix.  Search using the prefix instead.
            //
            // This fixes the case where we have these 4 inputs,
            // starting at "AB":
            //
            // 1: AB
            // 2: ABC
            // 3: ABB
            // 3: ABD
            //
            // If we search for "ABC", we'll get between "ABB" and "ABD" on
            // input 3.
            //
            // Instead, we notice that "AB" is a prefix, so we search again for
            // "AB" and get the correct position on input 3.

            target_seq = maybe_prefix;
            found_shorter = true;
            found_shorter_count++;
            NOTE_TEST_COVERAGE(make_mergemap);
            break;
          } else {
            pos = maybe_pos;
            NOTE_TEST_COVERAGE(make_mergemap);
          }
        }
      } else {
        NOTE_TEST_COVERAGE(make_mergemap);
      }
      indexes[input_id] = pos;
    }
  }
  NOTE_TEST_COVERAGE_IF(make_mergemap, found_shorter_count == 1);
  NOTE_TEST_COVERAGE_IF(make_mergemap, found_shorter_count == 2);
  NOTE_TEST_COVERAGE_IF(make_mergemap, found_shorter_count > 2);
}

void make_mergemap::add_to_queue(std::priority_queue<queue_entry>& queue,
                                 unsigned flat, size_t entry_id,
                                 size_t limit_idx) const {
  if (entry_id == limit_idx) {
    // No more entries to process from this part.
    return;
  }

  CHECK_LT(entry_id, limit_idx);
  CHECK_LT(entry_id, m_flats[flat]->size());

  queue_entry entry;
  entry.cur_slice = m_flats[flat]->get(entry_id);
  entry.flat = flat;
  entry.entry_id = entry_id;

  queue.push(entry);
}

void make_mergemap::count_range(size_t start, size_t limit) {
  dna_slice start_seq = seq_for_pos(start, empty_meaning::start);
  dna_slice limit_seq = seq_for_pos(limit, empty_meaning::limit);

  std::vector<size_t> start_idx(m_flats.size());
  std::vector<size_t> limit_idx(m_flats.size());
  std::vector<std::vector<bool>> bits(m_flats.size());
  std::priority_queue<queue_entry> queue;

  positions_for_seq(start_seq, start_idx, empty_meaning::start);
  positions_for_seq(limit_seq, limit_idx, empty_meaning::limit);

  for (unsigned i = 0; i < m_flats.size(); ++i) {
    add_to_queue(queue, i, start_idx[i], limit_idx[i]);
  }

  size_t merged_entries = 0;
  while (!queue.empty()) {
    // Extend merged_entries in blocks of 4k
    if ((merged_entries & 32767) == 0) {
      for (unsigned i = 0; i < m_flats.size(); ++i) {
        bits[i].resize(merged_entries + 32768, false /* initialize as unset */);
      }
    }

    // Extract the next entry in the merged seqset.
    queue_entry entry = queue.top();
    queue.pop();

    dna_slice slice = entry.cur_slice;
    bits[entry.flat][merged_entries] = 1;
    add_to_queue(queue, entry.flat, entry.entry_id + 1, limit_idx[entry.flat]);

    // Extract any entries that are duplicates of it
    while (!queue.empty() && is_equal_or_prefix(slice, queue.top().cur_slice)) {
      NOTE_TEST_COVERAGE_IF(make_mergemap,
                            slice.size() != queue.top().cur_slice.size());
      NOTE_TEST_COVERAGE_IF(make_mergemap,
                            slice.size() == queue.top().cur_slice.size());

      entry = queue.top();
      queue.pop();

      slice = entry.cur_slice;
      bits[entry.flat][merged_entries] = 1;
      add_to_queue(queue, entry.flat, entry.entry_id + 1,
                   limit_idx[entry.flat]);
    }

    merged_entries++;
  }

  chunk_result result;
  result.merged_entries = merged_entries;
  result.bits = std::move(bits);

  std::lock_guard<std::mutex> l(m_mu);
  CHECK(!m_chunk_results.count(start_seq));
  m_total_merged_entries += result.merged_entries;
  m_chunk_results[start_seq] = std::move(result);
}

void make_mergemap::fill_mergemap(unsigned input_id,
                                  seqset_mergemap_builder* mergemap,
                                  progress_handler_t progress) {
  SPLOG("Filling mergemap for part index %d", input_id);

  std::vector<std::pair<size_t, chunk_result*>> fill_positions;
  size_t index = 0;
  for (auto& r : m_chunk_results) {
    fill_positions.emplace_back(std::make_pair(index, &r.second));
    index += r.second.merged_entries;
  }
  CHECK_EQ(m_total_merged_entries, index);

  parallel_for(  //
      0, fill_positions.size(),
      [input_id, mergemap, &fill_positions](size_t chunk_idx) {
        size_t entry_index = fill_positions[chunk_idx].first;
        chunk_result* chunk = fill_positions[chunk_idx].second;
        std::vector<bool>* bits(&chunk->bits[input_id]);

        for (auto i = bits->begin(); i != bits->end(); ++i) {
          if (*i) {
            mergemap->set(entry_index);
          }
          entry_index++;
        }

        // Free up memory used.
        bits->clear();
      },
      progress);

  mergemap->finalize();
  SPLOG("Done bitcount for part index %d", input_id);
}

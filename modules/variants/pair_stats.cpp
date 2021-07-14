#include "modules/variants/pair_stats.h"

#include <random>

namespace variants {

void pair_stats::calc_stats(progress_handler_t progress) {
  std::vector<int64_t> offsets;

  std::mt19937_64 rand_source;
  CHECK_GT(m_seqset->size(), 0);
  std::uniform_int_distribution<uint64_t> rand_seqset_id(0, m_seqset->size() - 1);

  size_t cross_scaffold = 0;
  size_t tot_entries = 0;
  size_t loop_count = 0;
  size_t bad_direction = 0;

  while (offsets.size() < k_num_samples && loop_count < k_max_attempts) {
    ++loop_count;
    uint64_t orig_seqset_id = rand_seqset_id(rand_source);

    tot_entries++;

    auto ref_map_entry = m_ref_map->get(orig_seqset_id);
    if (!ref_map_entry.is_unique()) {
      continue;
    }

    auto read_id_range = m_readmap->entry_to_index(orig_seqset_id);
    if (read_id_range.first == read_id_range.second) {
      continue;
    }

    std::uniform_int_distribution<uint32_t> rand_read_id(read_id_range.first,
                                                         read_id_range.second - 1);
    uint32_t orig_read_id = rand_read_id(rand_source);

    if (!m_readmap->get_is_forward(orig_read_id)) {
      // TODO(nils): Remove this once readmaps actually store rev_comp
      // information.
      continue;
    }
    if (!m_readmap->has_mate(orig_read_id)) {
      continue;
    }

    uint32_t mate_read_id = m_readmap->get_mate(orig_read_id);
    if (!m_readmap->get_is_forward(mate_read_id)) {
      continue;
    }
    uint64_t mate_seqset_id = m_readmap->index_to_entry(mate_read_id);

    auto mate_ref_map_entry = m_ref_map->get(mate_seqset_id);
    if (!mate_ref_map_entry.is_unique()) {
      continue;
    }

    seq_position orig_loc;
    bool orig_loc_rc;
    get_ref_loc(orig_seqset_id, ref_map_entry, orig_loc, orig_loc_rc);

    seq_position mate_loc;
    bool mate_loc_rc;
    get_ref_loc(mate_seqset_id, mate_ref_map_entry, mate_loc, mate_loc_rc);

    if (orig_loc.scaffold_id != mate_loc.scaffold_id) {
      cross_scaffold++;
      continue;
    }

    if (orig_loc_rc == mate_loc_rc) {
      bad_direction++;
      continue;
    }

    int64_t distance = int64_t(mate_loc.position) - int64_t(orig_loc.position);
    if (orig_loc_rc) {
      distance = -distance;
    }
    offsets.push_back(distance);
    progress(offsets.size() * 1. / k_num_samples);
  }

  if (offsets.empty()) {
    SPLOG(
        "No pairs found to calculate distance statistics for (skipped %ld "
        "cross-extent and %ld mismatched direction pairs",
        cross_scaffold, bad_direction);
    m_median_pair_offset.reset();
    return;
  }

  auto median_it = offsets.begin() + (offsets.size() / 2);
  std::nth_element(offsets.begin(), median_it, offsets.end());
  m_median_pair_offset.emplace(*median_it);

  SPLOG(
      "After scanning %ld seqset entries, Found %ld pairs with a median "
      "distance of %ld (skipping %ld cross-extent pairs and %ld pairs with "
      "mismatched direction)",
      tot_entries, offsets.size(), *m_median_pair_offset, cross_scaffold, bad_direction);
}

void pair_stats::get_ref_loc(uint64_t seqset_id, const ref_map::entry& rmap_entry,
                             seq_position& loc, bool& loc_rc) {
  CHECK(rmap_entry.is_unique());
  if (rmap_entry.rev_match()) {
    loc_rc = true;
  } else {
    loc_rc = false;
  }

  dna_sequence seq = m_seqset->ctx_entry(seqset_id).sequence();
  if (loc_rc) {
    seq = seq.rev_comp();
  }

  bwt_range ref_range = m_ref->get_bwt().find(seq);
  CHECK_EQ(1, ref_range.matches()) << seq << " loc_rc: " << loc_rc;
  uint32_t flattened_pos = ref_range.get_match(0);

  loc = m_ref->get_seq_position(flattened_pos);
  if (loc_rc) {
    // Facing backwards, and we want the offset of the beginning of
    // the sequence.
    loc.position += seq.size();
  }
}

}  // namespace variants

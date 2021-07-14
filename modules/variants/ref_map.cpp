#include "modules/variants/ref_map.h"

#include "modules/io/parallel.h"

namespace variants {

constexpr unsigned ref_map::k_fwd_flag;
constexpr unsigned ref_map::k_rev_flag;
constexpr unsigned ref_map::k_count_mask;
constexpr size_t ref_map::k_min_chunk_size;
constexpr size_t ref_map::k_num_flush_buckets;
constexpr size_t ref_map::k_flush_bucket_size;

ref_map::ref_map(const seqset* the_seqset, const reference* ref)
    : m_seqset(the_seqset), m_ref(ref) {
  m_mutable_ref_map.reset(new mutable_packed_vector<unsigned, 8>(m_seqset->size(), "ref_map"));
  m_ref_map = m_mutable_ref_map;
}

ref_map::ref_map(const seqset* the_seqset, const reference* ref,
                 const spiral_file_create_state& state)
    : m_seqset(the_seqset), m_ref(ref) {
  m_mutable_ref_map.reset(new mutable_packed_vector<unsigned, 8>(state, m_seqset->size()));
  m_ref_map = m_mutable_ref_map;
}

ref_map::ref_map(const seqset* the_seqset, const reference* ref,
                 const spiral_file_open_state& state)
    : m_seqset(the_seqset), m_ref(ref) {
  m_ref_map.reset(new packed_vector<unsigned, 8>(state));
  CHECK_EQ(m_ref_map->size(), m_seqset->size());
}

namespace {

struct extent_slice {
  size_t prestart_len;
  dna_slice slice;
  bool is_rev_comp;
};

}  // namespace

void ref_map::flush_updates(std::vector<uint64_t>& seqset_ids, bool is_rev_comp) {
  for (size_t flush_id = 0; flush_id != k_num_flush_buckets; ++flush_id) {
    std::lock_guard<std::mutex> l(m_flush_bucket_mu[flush_id]);
    size_t flush_start = flush_id * m_seqset_entries_per_flush_bucket;
    size_t flush_limit = (flush_id + 1) * m_seqset_entries_per_flush_bucket;
    for (uint64_t seqset_id : seqset_ids) {
      if (seqset_id < flush_start || seqset_id >= flush_limit) {
        continue;
      }
      uint8_t old_val = m_mutable_ref_map->at(seqset_id);
      unsigned old_count = old_val & k_count_mask;
      unsigned new_count = old_count;
      if (old_count != k_count_mask) {
        new_count = old_count + 1;
      }
      uint8_t new_val = (old_val & ~k_count_mask) | new_count;

      if (is_rev_comp) {
        new_val |= k_fwd_flag;
      } else {
        new_val |= k_rev_flag;
      }
      m_mutable_ref_map->at(seqset_id).set_unlocked(new_val);
    }
  }
  seqset_ids.clear();
}

void ref_map::build(progress_handler_t progress) {
  const flat_ref& flat = m_ref->get_flat_ref();
  const flat_ref::index_t& index = flat.get_index();

  std::vector<extent_slice> ref_slices;
  size_t tot_ref_bases = 0;

  for (const flat_ref::extent_t& extent : index.extents) {
    tot_ref_bases += extent.size;
    auto data_start = m_ref->get_dna(extent.flat);
    dna_slice seq = dna_slice(data_start, data_start + extent.size);

    size_t chunk_size = seq.size() / 100;
    if (chunk_size < k_min_chunk_size) {
      chunk_size = k_min_chunk_size;
    }

    for (size_t chunk_start = 0; chunk_start < seq.size(); chunk_start += chunk_size) {
      size_t chunk_end = chunk_start + chunk_size;
      if (chunk_end > seq.size()) {
        chunk_end = seq.size();
      }

      size_t prestart = chunk_start;
      if (prestart >= 256) {
        prestart -= 256;
      } else {
        prestart = 0;
      }

      size_t prestart_rc = chunk_end;
      prestart_rc += 256;
      if (prestart_rc > seq.size()) {
        prestart_rc = seq.size();
      }

      extent_slice slice;
      slice.prestart_len = chunk_start - prestart;
      slice.slice = seq.subseq(prestart, chunk_end - prestart);
      slice.is_rev_comp = false;
      ref_slices.push_back(slice);

      extent_slice rc_slice;
      rc_slice.prestart_len = prestart_rc - chunk_end;
      rc_slice.slice = seq.subseq(chunk_start, prestart_rc - chunk_start).rev_comp();
      rc_slice.is_rev_comp = true;
      ref_slices.push_back(rc_slice);
    }
  }

  SPLOG(
      "Marking %ld bases in %ld extents (%ld extent sections, including RCs) "
      "as reference",
      tot_ref_bases, index.extents.size(), ref_slices.size());

  std::atomic<size_t> tot_marked{0};
  m_seqset_entries_per_flush_bucket =
      (m_seqset->size() + (k_num_flush_buckets - 1)) / k_num_flush_buckets;
  // Make sure it doesn't fall between bytes
  m_seqset_entries_per_flush_bucket += sizeof(uint64_t);
  m_seqset_entries_per_flush_bucket &= ~size_t(sizeof(uint64_t) - 1);
  CHECK_GE(m_seqset_entries_per_flush_bucket * k_num_flush_buckets, m_seqset->size());
  CHECK_EQ(m_seqset_entries_per_flush_bucket % sizeof(uint64_t), 0);
  parallel_for(  //
      0, ref_slices.size(),
      [this, &ref_slices, &tot_marked](uint64_t ref_slice_id) {
        size_t chunk_entries = 0;

        const extent_slice& slice = ref_slices[ref_slice_id];
        seqset_range r = m_seqset->ctx_begin();

        auto prestart_len_left = slice.prestart_len;

        auto pos = slice.slice.begin();
        auto limit = slice.slice.end();

        std::vector<uint64_t> seqset_ids;

        while (pos != limit) {
          r = r.push_front_drop((*pos).complement());
          CHECK(r.valid());
          ++pos;
          if (prestart_len_left) {
            --prestart_len_left;
            continue;
          }

          uint64_t seqset_id = r.begin();
          if (seqset_id + 1 != r.end()) {
            continue;
          }
          if (r.size() != m_seqset->entry_size(seqset_id)) {
            continue;
          }

          ++chunk_entries;
          seqset_ids.push_back(seqset_id);
          if (seqset_ids.size() >= k_flush_bucket_size) {
            flush_updates(seqset_ids, slice.is_rev_comp);
          }
        }
        flush_updates(seqset_ids, slice.is_rev_comp);
        CHECK(seqset_ids.empty());
        CHECK(!prestart_len_left);
        tot_marked.fetch_add(chunk_entries);
      },
      progress);
  SPLOG("%ld nodes marked by walking reference (%.2f%%)", tot_marked.load(),
        tot_marked.load() * 100. / m_seqset->size());
}

ref_map::entry ref_map::get(uint64_t seqset_id) const {
  CHECK_LT(seqset_id, m_ref_map->size());
  return entry(m_ref_map->at(seqset_id));
}

boost::optional<ref_anchor> ref_map::get_unique_ref_anchor(uint64_t seqset_id) const {
  auto rm_entry = get(seqset_id);

  if (!rm_entry.is_unique()) {
    return boost::none;
  }

  ref_anchor pos;

  if (rm_entry.rev_match()) {
    pos.rev_comp = true;
  } else {
    pos.rev_comp = false;
  }

  dna_sequence seq = m_seqset->ctx_entry(seqset_id).sequence();
  if (pos.rev_comp) {
    seq = seq.rev_comp();
  }

  bwt_range ref_range = m_ref->get_bwt().find(seq);
  CHECK_EQ(1, ref_range.matches()) << seq << " loc_rc: " << pos.rev_comp;
  uint32_t flattened_pos = ref_range.get_match(0);

  pos.pos = m_ref->get_seq_position(flattened_pos);
  if (pos.rev_comp) {
    // Facing backwards, and we want the offset of the beginning of
    // the sequence.
    pos.pos.position += seq.size();
  }
  return pos;
}

dna_slice ref_map::get_ref_slice(const ref_anchor& anchor) const {
  const auto& refasm = m_ref->get_assembly();
  size_t flat_pos = refasm.flatten(anchor.pos);
  const supercontig& sc = refasm.get_supercontig(flat_pos);

  if (anchor.rev_comp) {
    return dna_slice(m_ref->get_dna(sc.tot_offset), m_ref->get_dna(flat_pos)).rev_comp();
  } else {
    return dna_slice(m_ref->get_dna(flat_pos), m_ref->get_dna(sc.tot_offset + sc.len));
  }
}

}  // namespace variants

#include "modules/bio_base/refmatch.h"

#include "modules/io/parallel.h"

constexpr unsigned refmatch::k_fwd_flag;
constexpr unsigned refmatch::k_rev_flag;
constexpr unsigned refmatch::k_count_mask;

const product_version refmatch::k_refmatch_version{"1.0.0"};

refmatch::refmatch(const seqset* the_seqset, const reference* ref)
    : m_seqset(the_seqset), m_ref(ref) {}

refmatch::refmatch(const seqset* the_seqset, const reference* ref,
                   const spiral_file_open_state& state)
    : m_seqset(the_seqset), m_ref(ref) {
  state.enforce_max_version("refmatch", k_refmatch_version);
  m_per_entry.emplace(state.open_subpart("per-entry"));
  CHECK_EQ(m_per_entry->size(), m_seqset->size());
  membuf overflow_ids_membuf = state.open_membuf("overflow-ids");
  membuf overflow_counts_membuf = state.open_membuf("overflow-counts");
  unsigned n_overflow = overflow_ids_membuf.size() / sizeof(uint64_t);
  CHECK_EQ(n_overflow, overflow_counts_membuf.size() / sizeof(unsigned));
  const uint64_t* overflow_ids =
      reinterpret_cast<const uint64_t*>(overflow_ids_membuf.data());
  const unsigned* overflow_counts =
      reinterpret_cast<const unsigned*>(overflow_counts_membuf.data());

  m_overflow.reserve(n_overflow);

  for (size_t i = 0; i < n_overflow; ++i) {
    bool was_new_entry =
        m_overflow.insert(std::make_pair(overflow_ids[i], overflow_counts[i]))
            .second;
    CHECK(was_new_entry);
  }
}

refmatch::entry refmatch::get(uint64_t seqset_id) const {
  CHECK(m_per_entry);
  bool has_fwd = false;
  bool has_rev = false;

  uint8_t val = m_per_entry->at(seqset_id);

  if (val & k_fwd_flag) {
    has_fwd = true;
  }

  if (val & k_rev_flag) {
    has_rev = true;
  }

  unsigned count = val & k_count_mask;
  if (count == k_count_mask) {
    auto it = m_overflow.find(seqset_id);
    if (it != m_overflow.end()) {
      count += it->second;
    }
  }

  return entry(has_fwd, has_rev, count);
}

namespace {

struct extent_slice {
  size_t prestart_len;
  dna_slice slice;
  bool is_rev_comp;
};

}  // namespace

size_t refmatch_builder::g_min_chunk_size = 25600;

refmatch_builder::refmatch_builder(const seqset* the_seqset,
                                   const reference* ref)
    : refmatch(the_seqset, ref) {}

void refmatch_builder::build(const spiral_file_create_state& state,
                             progress_handler_t progress) {
  state.set_version("refmatch", k_refmatch_version);
  m_mutable_per_entry.emplace(state.create_subpart("per-entry"),
                              m_seqset->size());
  walk_reference(progress);
  save_overflow(state);
}

void refmatch_builder::walk_reference(progress_handler_t progress) {
  const flat_ref& flat = m_ref->get_flat_ref();
  const flat_ref::index_t& index = flat.get_index();

  std::vector<extent_slice> ref_slices;
  size_t tot_ref_bases = 0;

  for (const flat_ref::extent_t& extent : index.extents) {
    tot_ref_bases += extent.size;
    auto data_start = m_ref->get_dna(extent.flat);
    dna_slice seq = dna_slice(data_start, data_start + extent.size);

    size_t chunk_size = seq.size() / 100;
    if (chunk_size < g_min_chunk_size) {
      chunk_size = g_min_chunk_size;
    }

    for (size_t chunk_start = 0; chunk_start < seq.size();
         chunk_start += chunk_size) {
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
      rc_slice.slice =
          seq.subseq(chunk_start, prestart_rc - chunk_start).rev_comp();
      rc_slice.is_rev_comp = true;
      ref_slices.push_back(rc_slice);
    }
  }

  SPLOG(
      "Marking %ld bases in %ld extents (%ld extent sections, including RCs) "
      "as reference",
      tot_ref_bases, index.extents.size(), ref_slices.size());

  size_t tot_marked = 0;
  std::mutex mu;
  parallel_for(  //
      0, ref_slices.size(),
      [this, &ref_slices, &tot_marked, &mu](uint64_t ref_slice_id) {
        size_t chunk_entries = 0;
        std::unordered_map<uint64_t, unsigned> chunk_overflow;

        const extent_slice& slice = ref_slices[ref_slice_id];
        seqset_range r = m_seqset->ctx_begin();

        auto prestart_len_left = slice.prestart_len;

        auto pos = slice.slice.begin();
        auto limit = slice.slice.end();

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

          for (;;) {
            uint8_t old_val = m_mutable_per_entry->at(seqset_id);
            unsigned old_count = old_val & k_count_mask;
            unsigned new_count = old_count;
            bool did_overflow = false;
            if (old_count == k_count_mask) {
              did_overflow = true;
            } else {
              new_count = old_count + 1;
            }
            uint8_t new_val = (old_val & ~k_count_mask) | new_count;

            if (slice.is_rev_comp) {
              new_val |= k_fwd_flag;
            } else {
              new_val |= k_rev_flag;
            }

            if (m_mutable_per_entry->at(seqset_id).compare_and_swap(old_val,
                                                                    new_val)) {
              ++chunk_entries;
              if (did_overflow) {
                chunk_overflow[seqset_id]++;
              }
              break;
            }
          }
        }

        CHECK(!prestart_len_left);

        std::lock_guard<std::mutex> l(mu);
        tot_marked += chunk_entries;
        for (const auto& o : chunk_overflow) {
          m_overflow[o.first] += o.second;
        }
      },
      progress);
  SPLOG("%ld nodes marked (%.2f%%) including %ld overflow entries", tot_marked,
        tot_marked * 100. / m_seqset->size(), m_overflow.size());
}

void refmatch_builder::save_overflow(const spiral_file_create_state& state) {
  size_t n_overflow = m_overflow.size();

  mutable_membuf overflow_ids_membuf =
      state.create_membuf("overflow-ids", sizeof(uint64_t) * n_overflow);
  mutable_membuf overflow_counts_membuf =
      state.create_membuf("overflow-counts", sizeof(unsigned) * n_overflow);
  uint64_t* overflow_ids =
      reinterpret_cast<uint64_t*>(overflow_ids_membuf.mutable_data());
  unsigned* overflow_counts =
      reinterpret_cast<unsigned*>(overflow_counts_membuf.mutable_data());

  auto it = m_overflow.begin();
  for (size_t i = 0; i < n_overflow; ++i) {
    CHECK(it != m_overflow.end());
    overflow_ids[i] = it->first;
    overflow_counts[i] = it->second;

    ++it;
  }
  CHECK(it == m_overflow.end());
}

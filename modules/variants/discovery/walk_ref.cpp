#include "modules/variants/discovery/walk_ref.h"

#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/variants/discovery/path.h"
#include "modules/variants/discovery/push_search.h"

namespace variants {
namespace discovery {

static constexpr bool k_dbg = false;

void walk_ref_t::walk_ref(aoffset_t start, aoffset_t limit) {
  for (const auto& ext : m_view->get_scaffold().extents()) {
    if (ext.offset >= limit) {
      break;
    }

    aoffset_t cur_offset = ext.offset + ext.sequence.size();
    if (cur_offset < start) {
      continue;
    }

    auto it = ext.sequence.rcbegin();
    if (cur_offset > aoffset_t(limit + opts().readmap->max_read_len())) {
      aoffset_t diff = cur_offset - aoffset_t(limit + opts().readmap->max_read_len());
      CHECK_GT(diff, 0);
      it += diff;
      cur_offset -= diff;
    }

    // Divide by 2 to avoid overflow when doing arithmetic.
    aoffset_t bases_since_read = std::numeric_limits<aoffset_t>::max() / 2;
    seqset_range r = opts().seqset->ctx_begin();
    wr_range_info_t* prev_wri = nullptr;
    while (it != ext.sequence.rcend()) {
      dna_base ref_base = (*it).complement();
      seqset_range new_r = r.push_front_drop(ref_base);
      if (new_r.size() != r.size() + 1) {
        // Dropped some bases; save the ref location.
        if (prev_wri) {
          save_ref_range(prev_wri);
        }
        prev_wri = nullptr;
      }
      r = new_r;
      --cur_offset;
      ++it;

      if (opts().bidir_treat_ref_as_reads) {
        bases_since_read = 0;
      } else {
        boost::optional<uint32_t> longest = opts().readmap->get_longest_prefix_read_id(r);
        if (longest) {
          bases_since_read = 0;
        } else {
          ++bases_since_read;
        }
      }

      if (r.size() < opts().min_overlap) {
        continue;
      }

      if (cur_offset >= limit) {
        continue;
      }

      if (cur_offset < start) {
        break;
      }

      boost::optional<dna_base> next_ref_base{};
      if (it != ext.sequence.rcend()) {
        next_ref_base = (*it).complement();
      }
      prev_wri = add_ref_range(cur_offset, r, dna_slice(it - r.size(), it).rev_comp(),
                               next_ref_base, bases_since_read);
    }
    if (prev_wri) {
      save_ref_range(prev_wri);
      prev_wri = nullptr;
    }
  }
}

void walk_ref_t::save_ref_range(wr_range_info_t* wri) {
  range_info_t& ri = m_view->range_info()[wri->r];
  ri.reference_offsets.push_back(wri->offset);
  wri->saved = true;
}

walk_ref_t::wr_range_info_t* walk_ref_t::add_ref_range(aoffset_t offset, const seqset_range& r,
                                                       dna_slice seq,
                                                       boost::optional<dna_base> next_ref_base,
                                                       int bases_since_read) {
  m_seen_ranges.emplace_back(r);
  wr_range_info_t wri;
  wri.r = r;
  wri.seq = seq;
  wri.next_ref_base = next_ref_base;
  wri.offset = offset;
  wri.bases_since_read = bases_since_read;
  if (opts().bidir_validate_trace_state > 1) {
    CHECK_EQ(wri.r.sequence(), wri.seq);
    CHECK_EQ(m_view->get_scaffold().subscaffold(offset, r.size()).get_simple(), seq);
  }
  m_wr_range_info.emplace_back(std::move(wri));
  return &m_wr_range_info.back();
}

void walk_ref_t::check_invariants() const {
  for (const auto& wri : m_wr_range_info) {
    const auto& ri = m_view->range_info().at(wri.r);

    if (m_view->opts().bidir_validate_trace_state > 1) {
      CHECK_EQ(wri.r.sequence(), wri.seq);
    }
    CHECK_EQ(wri.seq, m_view->get_scaffold().subscaffold(wri.offset, wri.r.size()).get_simple());
    auto offsets = ri.reference_offsets;
    if (wri.saved) {
      CHECK(std::find(offsets.begin(), offsets.end(), wri.offset) != offsets.end())
          << "Offset " << wri.offset << " not present for range " << wri.seq << "?";
    }
  }
}

void walk_ref_t::init_pairs_and_push() {
  std::sort(m_seen_ranges.begin(), m_seen_ranges.end(), seqset_range_comparer());
  if (k_dbg) {
    seqset_range prev_r;
    unsigned unique_seen = 0;
    for (const auto& r : m_seen_ranges) {
      if (r != prev_r) {
        prev_r = r;
        ++unique_seen;
      }
    }
    std::cout << m_seen_ranges.size() << " seen counts entries, " << unique_seen << " unique\n";
  }
  for (const wr_range_info_t& wri : m_wr_range_info) {
    bool ref_is_unique = true;
    bool matched_our_range = false;
    auto it = std::lower_bound(m_seen_ranges.begin(), m_seen_ranges.end(), wri.r,
                               seqset_range_comparer());
    while (it != m_seen_ranges.end()) {
      if (it->end() > wri.r.end()) {
        break;
      }
      if (matched_our_range) {
        // Found another match besides this range info entry.
        ref_is_unique = false;
        break;
      }
      matched_our_range = true;
      ++it;
    }
    CHECK(matched_our_range);

    if (ref_is_unique || opts().trace_ambiguous_ref) {
      if (aoffset_t(wri.r.size()) >= aoffset_t(m_view->opts().min_overlap) + wri.bases_since_read) {
        path p(opts().readmap, wri.seq, wri.r, wri.r.size() /* path overlap */,
               wri.bases_since_read, wri.r.size() /* anchor length */);
        m_view->add_push_traces(p, wri.offset, wri.next_ref_base);
      }
    }

    if (ref_is_unique || !opts().ignore_ambiguous_ref_pairs) {
      m_view->add_pair_offset_support_for_range(wri.offset, wri.offset, wri.r);
    }
  }

  m_wr_range_info.clear();
  m_seen_ranges.clear();
}

}  // namespace discovery
}  // namespace variants

#include "modules/variants/discovery/view.h"
#include "modules/bio_base/readmap.h"
#include "modules/variants/discovery/branch.h"
#include "modules/variants/discovery/pop_search.h"
#include "modules/variants/discovery/push_search.h"
#include "modules/variants/discovery/rejoin.h"

static constexpr bool k_dbg = false;
namespace variants {
namespace discovery {

right_partial::right_partial(dna_slice new_seq, aoffset_t new_outer_right_offset,
                             unsigned new_pair_match_count)
    : seq(new_seq),
      outer_right_offset(new_outer_right_offset),
      pair_match_count(new_pair_match_count) {}

std::pair<std::unique_ptr<view_t>, std::unique_ptr<view_t>> view_t::create_view_pair(
    const scaffold& s, state* st) {
  std::pair<std::unique_ptr<view_t>, std::unique_ptr<view_t>> result;
  auto& fwd = result.first;
  auto& rev = result.second;

  fwd.reset(new view_t);
  rev.reset(new view_t);

  fwd->m_state = st;
  rev->m_state = st;

  fwd->m_is_rev_comp = false;
  rev->m_is_rev_comp = true;

  fwd->m_reverse = rev.get();
  rev->m_reverse = fwd.get();

  fwd->m_scaffold = s;
  rev->m_scaffold = s.rev_comp();

  return result;
}

offset_info view_t::get_offset_info(aoffset_t offset, bool fwd) const {
  offset_info result;
  if (m_is_rev_comp) {
    result = m_state->get_offset_info(reverse_offset(offset), !fwd);
    if (result.ploids_remaining > 0) {
      result.ref_remaining_limit = reverse_offset(result.ref_remaining_limit);
    }
  } else {
    result = m_state->get_offset_info(offset, fwd);
  }
  if (result.ploids_remaining > 0) {
    if (fwd) {
      CHECK_GE(result.ref_remaining_limit, offset);
    } else {
      CHECK_LE(result.ref_remaining_limit, offset);
    }
    CHECK_GE(result.ref_remaining_limit, 0);
    CHECK_LE(result.ref_remaining_limit, m_scaffold.end_pos());
  }
  return result;
}

void view_t::add_right_partial(const seqset_range& r, right_partial rp) {
  m_range_info[r].right_partials.emplace_back(std::move(rp));
}

void view_t::check_invariants() const {
  for (const auto& r_and_ri : m_range_info) {
    const auto& r = r_and_ri.first;
    const auto& ri = r_and_ri.second;

    for (aoffset_t ref_offset : ri.reference_offsets) {
      auto ext = m_scaffold.split_extent_at(ref_offset);
      CHECK_GE(ext.second.size(), r.size());
      CHECK_EQ(ext.second.subseq(0, r.size()), r.sequence());
    }

    for (const auto& rp : ri.right_partials) {
      unsigned shared = shared_ref_bases_to_left(rp.outer_right_offset, rp.seq);
      CHECK_GE(shared, opts().bidir_min_anchor_len) << rp;
    }
  }

  for (const auto& br_item : m_branches) {
    const auto& br_key = br_item.first;
    const auto& br = br_item.second;
    br->check_invariants();
    CHECK_EQ(br->right_push_view_offset(), br_key.second);
    CHECK_EQ(br->first_base(), br_key.first);
    std::string bstr;
    bstr += char(br->first_base());
    CHECK_NE(get_scaffold().subscaffold_str(br->right_push_view_offset() - 1, 1), bstr);
  }
}

unsigned view_t::shared_ref_bases_to_right(aoffset_t ref_offset, dna_slice slice) const {
  auto ext = get_scaffold().split_extent_at(ref_offset);
  return ext.second.shared_prefix_length(slice);
}

unsigned view_t::shared_ref_bases_to_left(aoffset_t ref_offset, dna_slice slice) const {
  return reverse_view()->shared_ref_bases_to_right(reverse_offset(ref_offset), slice.rev_comp());
}

void view_t::add_pair_offset_support(aoffset_t min_offset, aoffset_t max_offset,
                                     const readmap::read& rd) {
  CHECK_GE(max_offset, min_offset);

  if (!rd.has_mate()) {
    return;
  }

  seqset_range mate_r = rd.get_mate().get_rev_comp().get_seqset_entry();

  interval_t pair_valid;
  bool rd_faces_inward = opts().forward_pairs_face_inward ? rd.is_original_orientation()
                                                          : !rd.is_original_orientation();

  if (rd_faces_inward) {
    min_offset += opts().min_pair_distance;
    max_offset += opts().max_pair_distance;
  } else {
    min_offset -= opts().min_pair_distance;
    max_offset -= opts().max_pair_distance;
    if (min_offset > max_offset) {
      std::swap(min_offset, max_offset);
    }
  }

  if (rd_faces_inward) {
    // Store acceptable ranges for the beginning of mate_r, not the end.
    min_offset -= mate_r.size();
    max_offset -= mate_r.size();
  } else {
    // Min/Max distances should be measured from the end of rd, not the beginning.
    min_offset += rd.size();
    max_offset += rd.size();
  }

  m_range_info[mate_r].pair_supported_offsets += interval_t(min_offset, max_offset);
}

void view_t::add_pair_offset_support_for_range(aoffset_t min_offset, aoffset_t max_offset,
                                               const seqset_range& r) {
  CHECK_GE(max_offset, min_offset);

  size_t num_reads = 0;
  for (const readmap::read& rd __attribute__((unused)) : opts().readmap->get_prefix_reads(r)) {
    ++num_reads;
    if (num_reads > opts().max_pairs_per_read) {
      // Too many reads; don't populate pairing data
      return;
    }
  }

  CHECK_LE(num_reads, opts().max_pairs_per_read);

  for (const readmap::read& rd : opts().readmap->get_prefix_reads(r)) {
    add_pair_offset_support(min_offset, max_offset, rd);
  }
}

void view_t::add_push_traces(const path& p, aoffset_t right_offset,
                             boost::optional<dna_base> base_to_skip, branch* on_branch) {
  if (!on_branch) {
    CHECK_EQ(p.anchor_len(), p.seq().size());
  } else {
    CHECK_GT(p.seq().size(), p.anchor_len());
    unsigned shared = shared_ref_bases_to_left(right_offset + p.anchor_len(), p.seq());
    CHECK_EQ(shared, p.anchor_len()) << "\nPath: " << p << "\nRight_offset: " << right_offset
                                     << "\n";
  }

  int size_needed = aoffset_t(opts().min_overlap) + p.bases_since_read();
  CHECK_GT(size_needed, 0);
  if (int(p.range().size()) < size_needed) {
    return;
  }
  for (dna_base b : dna_bases()) {
    if (base_to_skip && b == *base_to_skip) {
      // Not a branch
      continue;
    }

    seqset_range pushed_r = p.range().push_front_drop(b, size_needed);
    if (!pushed_r.valid()) {
      continue;
    }
    path pushed = p;
    pushed.push_front_drop(b, pushed_r);

    branch* br = on_branch;
    if (!br) {
      br = get_branch(b, right_offset);
    }

    std::unique_ptr<push_search_entry> push =
        make_unique<push_search_entry>(std::move(pushed), 0 /* no pair matches yet */);
    if (k_dbg) {
      std::cout << "Generated push trace: " << push->describe(br) << "\n";
    }
    br->add_search_entry(std::move(push));
  }
}

void view_t::walk_assembly_variants(unsigned path_overlap, aoffset_t left_offset,
                                    aoffset_t left_anchor_len, aoffset_t right_offset,
                                    aoffset_t right_anchor_len, dna_slice seq, branch* on_branch) {
  if (k_dbg) {
    std::cout << "\nWalking assembly variants for " << seq
              << " left anchor len = " << left_anchor_len
              << " right anchor len = " << right_anchor_len << "\n";
  }
  CHECK_GE(right_offset + right_anchor_len, left_offset - left_anchor_len);

  CHECK_LT(right_anchor_len, seq.size());
  CHECK_GT(right_anchor_len, opts().bidir_min_anchor_len);
  CHECK_LT(left_anchor_len, seq.size());
  CHECK_GT(left_anchor_len, opts().bidir_min_anchor_len);

  if (opts().bidir_validate_trace_state) {
    CHECK_EQ(shared_ref_bases_to_right(left_offset - left_anchor_len, seq), left_anchor_len);
    CHECK_EQ(shared_ref_bases_to_left(right_offset + right_anchor_len, seq), right_anchor_len);
  }
  auto rc_it = seq.rcbegin();
  aoffset_t seq_bases_left = seq.size();

  // Traverse the right anchor to find the start of the variant
  seqset_range init_r = opts().seqset->ctx_begin();
  for (int i = 0; i < right_anchor_len; ++i) {
    CHECK(seq_bases_left);
    CHECK(rc_it != seq.rcend());

    init_r = init_r.push_front_drop((*rc_it).complement());
    --seq_bases_left;
    ++rc_it;
  }

  dna_slice right_anchor = dna_slice(seq).subseq(seq.size() - right_anchor_len, right_anchor_len);
  path p(opts().readmap, right_anchor, init_r, right_anchor_len, 0, right_anchor_len);

  // Add traces from the middle of the variant.
  while (seq_bases_left > left_anchor_len) {
    CHECK(seq_bases_left);
    CHECK(rc_it != seq.rcend());

    p.push_front_drop((*rc_it).complement());

    ++rc_it;
    --seq_bases_left;

    CHECK(rc_it != seq.rcend());
    CHECK(seq_bases_left);

    boost::optional<dna_base> next_b;
    CHECK(rc_it != seq.rcend());
    next_b.emplace((*rc_it).complement());

    add_push_traces(p, right_offset, next_b, on_branch);
    add_pair_offset_support_for_range(left_offset - left_anchor_len,
                                      right_offset + right_anchor_len, p.range());
  }

  CHECK_GT(seq_bases_left, 0);
  CHECK_LE(seq_bases_left, left_anchor_len);

  auto ext = get_scaffold().split_extent_at(left_offset - left_anchor_len + seq_bases_left);
  dna_slice rc_ref_slice_left = ext.first.rev_comp();
  auto rc_ref_it = rc_ref_slice_left.begin();
  auto rc_ref_it_end = rc_ref_slice_left.end();

  dna_slice bases_left_ref = rc_ref_slice_left.subseq(0, seq_bases_left).rev_comp();
  dna_slice bases_left = seq.subseq(0, seq_bases_left);
  CHECK_EQ(bases_left_ref, bases_left);

  CHECK(rc_ref_it == bases_left_ref.rcbegin());
  CHECK(rc_it == bases_left.rcbegin());

  unsigned ref_bases_pushed = 0;
  if (seq_bases_left < left_anchor_len) {
    ref_bases_pushed = left_anchor_len - seq_bases_left;
  }

  while (rc_ref_it != rc_ref_it_end) {
    p.push_front_drop((*rc_ref_it).complement());
    ++ref_bases_pushed;
    ++rc_ref_it;

    if (p.range().size() <= ref_bases_pushed) {
      break;
    }

    boost::optional<dna_base> next_b;
    if (rc_ref_it != rc_ref_it_end) {
      next_b.emplace((*rc_ref_it).complement());
    }

    add_push_traces(p, right_offset, next_b, on_branch);
    add_pair_offset_support_for_range(left_offset - left_anchor_len,
                                      right_offset + right_anchor_len, p.range());
  }

  if (k_dbg) {
    std::cout << "After walking variant for push traces in left anchor part of sequence " << seq
              << ", final path is: " << p << "\n";
  }
}

std::ostream& operator<<(std::ostream& os, const range_info_t& ri) {
  os << "\nRange info:\n";
  os << "  Reference offsets: ";
  for (aoffset_t offset : ri.reference_offsets) {
    os << " " << offset << "\n";
  }
  os << "  Pair supported offsets: " << ri.pair_supported_offsets << "\n";
  os << "  Partials:\n";
  for (const auto& rp : ri.right_partials) {
    os << "     " << rp << "\n";
  }
  return os;
}

bool view_t::has_ploids_remaining(aoffset_t left_offset, aoffset_t right_offset) const {
  offset_info oi = get_offset_info(left_offset, true /* fwd */);
  if (oi.ploids_remaining <= 0 || oi.ref_remaining_limit < right_offset) {
    return false;
  }
  return true;
}

branch* view_t::get_branch(dna_base b, aoffset_t right_offset) {
  if (opts().bidir_validate_trace_state) {
    auto ext = get_scaffold().split_extent_at(right_offset - 1);
    if (!ext.second.empty()) {
      CHECK_NE(b, ext.second[0]) << " branch at " << right_offset << " of "
                                 << (is_rev_comp() ? "fwd" : "rev") << " view of scaffold "
                                 << opts().scaffold_name;
    }
  }

  auto& branch_ptr = m_branches[branch_key(b, right_offset)];
  if (!branch_ptr) {
    constexpr bool k_show_new_branches = false;
    if (k_show_new_branches) {
      if (is_rev_comp()) {
        std::cout << "new rev branch at " << reverse_offset(right_offset) << " -> "
                  << b.complement() << "\n";
      } else {
        std::cout << "new fwd branch at " << b << " <- " << right_offset << "\n";
      }
    }
    branch_ptr = std::make_shared<branch>(this, b, right_offset);
    if (opts().bidir_validate_trace_state) {
      branch_ptr->check_invariants();
    }
  }
  return branch_ptr.get();
}

void view_t::discard_search_entries() {
  for (const auto& item : m_branches) {
    item.second->clear();
  }
}

std::vector<branch*> view_t::branches() const {
  std::vector<branch*> result;
  for (const auto& brinfo : m_branches) {
    result.push_back(brinfo.second.get());
  }
  return result;
}

}  // namespace discovery
}  // namespace variants

#include "modules/variants/discovery/rejoin.h"

namespace variants {
namespace discovery {

static constexpr bool k_trace = false;

const char rejoin_search_entry::k_tracer_name[] = "BIDIR";

rejoin_search_entry::rejoin_search_entry(unsigned path_overlap, aoffset_t left_offset,
                                         aoffset_t left_anchor_len, path p,
                                         unsigned pair_match_count)
    : branch_search_entry(search_entry_key(search_priority::REJOIN, p, pair_match_count)),
      m_left_offset(left_offset),
      m_left_anchor_len(left_anchor_len),
      m_path(std::move(p)) {
  if (m_key.path_overlap > m_left_anchor_len) {
    m_key.path_overlap = m_left_anchor_len;
  }

  if (m_path.range().size() > m_left_anchor_len) {
    unsigned drop_base_count = m_path.range().size() - m_left_anchor_len;
    m_key.tot_overlap_bases += drop_base_count;
    m_key.tot_overlap += m_left_anchor_len * drop_base_count;
  }
}

void rejoin_search_entry::notify_discard(branch* br) {
  if (br->push_view()->opts().bidir_tracer_emit_all_rejoins) {
    output(br, false /* don't generate more push traces if we're only
                      * emitting due to
                      * bidir_tracer_emit_all_rejoins */);
  }
}

search_result rejoin_search_entry::search_internal(branch* br) {
  view_t* v = br->push_view();
  bool left_within_range = m_left_offset >= br->push_view_farthest_left_offset();

  if (!v->opts().bidir_tracer_emit_all_rejoins && !left_within_range) {
    return search_result::STOP_SEARCHING;
  }

  output(br, left_within_range /* generate more push traces if we're
                                * within range.  */);
  return search_result::STOP_SEARCHING;
}

void rejoin_search_entry::output(branch* br, bool walk_more_vars) {
  if (k_trace || br->trace_enabled(m_path)) {
    std::cout << "Outputting rejoin " << describe(br) << "(walk_more=" << walk_more_vars << ")\n";
  }

  view_t* v = br->push_view();

  assembly_ptr a = make_unique<assembly>();
  a->tags.insert(k_tracer_name);
  a->assembly_id = allocate_assembly_id();
  a->min_overlap = m_key.path_overlap;

  aoffset_t left_offset = m_left_offset;
  aoffset_t right_offset = br->right_push_view_offset();

  dna_slice seq = m_path.seq();
  dna_slice rc_seq = seq.rev_comp();

  aoffset_t left_anchor_len = m_left_anchor_len;
  aoffset_t right_anchor_len = m_path.anchor_len();
  aoffset_t tot_anchor_len = left_anchor_len + right_anchor_len;

  aoffset_t anchor_overlap =
      std::max<aoffset_t>(tot_anchor_len - aoffset_t(seq.size()), left_offset - right_offset);

  aoffset_t rc_right_offset = v->reverse_offset(left_offset);
  aoffset_t rc_left_offset = v->reverse_offset(right_offset);

  if (v->is_rev_comp()) {
    a->left_offset = rc_left_offset - right_anchor_len;
    a->right_offset = rc_right_offset + left_anchor_len;
    a->left_anchor_len = right_anchor_len;
    a->right_anchor_len = left_anchor_len;
    a->seq = rc_seq;
  } else {
    a->left_offset = left_offset - left_anchor_len;
    a->right_offset = right_offset + right_anchor_len;
    a->left_anchor_len = left_anchor_len;
    a->right_anchor_len = right_anchor_len;
    a->seq = seq;
  }

  if (anchor_overlap > 0) {
    a->left_anchor_len -= anchor_overlap;
    a->right_anchor_len -= anchor_overlap;
    if (a->left_anchor_len <= 0 || a->right_anchor_len <= 0) {
      // Overlaps too much.
      return;
    }
  }

  check_assembly(*a, "discovery:rejoin");

  br->note_output(seq);
  br->get_state()->output_assembly(std::move(a), walk_more_vars);

  if (walk_more_vars) {
    view_t* rc_v = v->reverse_view();
    dna_base rc_branch_base = seq[m_left_anchor_len].complement();
    branch* rc_br = rc_v->get_branch(rc_branch_base, v->reverse_offset(m_left_offset));

    rc_br->note_output(seq.rev_comp());

    if (br->any_trace_enabled() || rc_br->any_trace_enabled()) {
      std::cout << "Clearing branches due to rejoin: " << describe(br) << "\n";
    }
    // Stop searching for the path that we just output.
    br->notify_rejoin(br, this);
    rc_br->notify_rejoin(br, this);

    // But continue searching for further variants that are different
    // than this path, if any exist.
    v->walk_assembly_variants(m_key.path_overlap, m_left_offset, m_left_anchor_len, right_offset,
                              right_anchor_len, seq, br);
    rc_v->walk_assembly_variants(m_key.path_overlap, rc_left_offset, right_anchor_len,
                                 rc_right_offset, m_left_anchor_len, rc_seq, rc_br);
  }
}

void rejoin_search_entry::check_invariants(const branch* br) const {
  CHECK_LE(m_key.path_overlap, m_path.path_overlap()) << describe(br);
  if (br->opts().bidir_validate_trace_state > 1) {
    br->check_path_invariants(m_path);
  }
  view_t* v = br->push_view();
  aoffset_t right_offset = br->right_push_view_offset();
  aoffset_t right_anchor_len = m_path.anchor_len();

  CHECK_LE(m_left_offset - m_left_anchor_len, right_offset + m_path.path_overlap()) << describe(br);

  CHECK_LT(m_left_anchor_len, m_path.size()) << describe(br);
  CHECK_LT(right_anchor_len, m_path.size()) << describe(br);

  unsigned shared_left =
      v->shared_ref_bases_to_right(m_left_offset - m_left_anchor_len, m_path.seq());
  CHECK_EQ(shared_left, m_left_anchor_len)
      << describe(br) << "\nRef before:\n"
      << v->get_scaffold().subscaffold_str(m_left_offset - m_left_anchor_len, m_left_anchor_len)
      << "\nRef after: " << v->get_scaffold().subscaffold_str(m_left_offset, 300)
      << "\nSeq: " << m_path.seq() << "\nLeft anchor len: " << m_left_anchor_len;

  unsigned shared_right =
      v->shared_ref_bases_to_left(right_offset + right_anchor_len, m_path.seq());
  CHECK_EQ(shared_right, right_anchor_len)
      << describe(br) << "\nRef before:\n"
      << v->get_scaffold().subscaffold_str(right_offset - 100, 100)
      << "\nRef after: " << v->get_scaffold().subscaffold_str(right_offset, 300)
      << "\nSeq: " << m_path.seq() << "\nRight anchor len: " << right_anchor_len;
}

std::string rejoin_search_entry::describe_internal(const branch* br) const {
  aoffset_t right_offset = br->right_push_view_offset();
  aoffset_t right_anchor_len = m_path.anchor_len();
  std::stringstream result;

  aoffset_t seq_len =
      aoffset_t(m_path.seq().size()) - right_anchor_len - aoffset_t(m_left_anchor_len);
  aoffset_t ref_len = right_offset - m_left_offset;
  aoffset_t svlen = seq_len - ref_len;

  view_t* v = br->push_view();

  if (v->is_rev_comp()) {
    result << "rev-rejoin@" << v->opts().scaffold_name << ":" << v->reverse_offset(right_offset)
           << "(al=" << right_anchor_len << ")->" << v->reverse_offset(m_left_offset)
           << "(al=" << m_left_anchor_len << ") svlen=" << svlen << " path=" << m_path;
  } else {
    result << "fwd-rejoin@" << v->opts().scaffold_name << ",ol=" << m_key.path_overlap << ":"
           << right_offset << "(al=" << right_anchor_len << ")->" << m_left_offset
           << "(al=" << m_left_anchor_len << ") svlen=" << svlen << " path=" << m_path;
  }

  return result.str();
}

}  // namespace discovery
}  // namespace variants

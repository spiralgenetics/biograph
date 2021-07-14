#include "modules/variants/discovery/branch.h"
#include "modules/variants/discovery/path.h"
#include "modules/variants/discovery/rejoin.h"
#include "modules/variants/discovery/view.h"

namespace variants {
namespace discovery {

static constexpr bool k_dbg = false;

search_result branch_search_entry::search(branch* br) { return search_internal(br); }

struct branch_search_entry_comparer {
  bool operator()(const branch_search_entry_ptr& lhs, const branch_search_entry_ptr& rhs) const {
    // Returns true if rhs is better than lhs.
    // Exhaust all push and rejoin searches first before trying POP searches.
    const auto& akey = lhs->get_key();
    const auto& bkey = rhs->get_key();
    return (*this)(akey, bkey);
  }

  bool operator()(const search_entry_key& akey, const search_entry_key& bkey) const {
    return akey < bkey;
  }
};

void branch::check_invariants() const {
  unsigned tracable_count = 0;
  if (!m_trace.empty()) {
    for (const auto& e : m_search_entries) {
      if (trace_enabled_for_entry(e.get())) {
        ++tracable_count;
      }
      e->check_invariants(this);
    }
  }

  CHECK_EQ(tracable_count, m_tracable_entry_count);

  // Make sure we really do split off of reference, and in a valid offset.
  auto splits = push_view()->get_scaffold().split_extent_at(m_right_offset);
  auto left_of_anchor = splits.first;
  if (!left_of_anchor.empty()) {
    CHECK_NE(m_first_base, left_of_anchor[left_of_anchor.size() - 1]);
  }
}

std::string branch::describe() const {
  std::stringstream out;
  if (m_push_view->is_rev_comp()) {
    out << "Branch " << opts().scaffold_name << ":" << left_pop_view_offset()
        << "(rev=" << right_push_view_offset() << ") -> " << m_first_base.complement();
  } else {
    out << "Branch " << m_first_base << " -> " << opts().scaffold_name << ":"
        << right_push_view_offset() << "(rev=" << left_pop_view_offset() << ")";
  }
  out << ", " << m_search_entries.size() << " entries\n";

  constexpr bool k_display_all_search_entries = false;
  if (k_display_all_search_entries) {
    for (const auto& e : m_search_entries) {
      out << "  " << e->describe(this) << "\n";
    }
  }
  return out.str();
}

view_t* branch::push_view() const { return m_push_view; }
view_t* branch::pop_view() const { return m_push_view->reverse_view(); }
state* branch::get_state() const { return m_push_view->get_state(); }

branch::branch(view_t* push_view, dna_base first_base, aoffset_t right_offset)
    : m_push_view(push_view), m_first_base(first_base), m_right_offset(right_offset) {
  m_steps_left = opts().bidir_max_branch_steps;
}

bool branch::update_ploids_remaining() {
  offset_info oi = push_view()->get_offset_info(m_right_offset, false /* !fwd */);
  if (oi.ploids_remaining <= 0) {
    m_ref_remaining = std::numeric_limits<aoffset_t>::max();
    return false;
  }
  m_ref_remaining = oi.ref_remaining_limit;
  return true;
}

void branch::notify_rejoin(branch* other_br, branch_search_entry* other_e) {
  update_ploids_remaining();
  // If we got a rejoin, discard all our intermediate searches.
  if (m_tracable_entry_count) {
    std::cout << "Clearing branch " << *this << " because of " << other_e->describe(other_br)
              << "\n";
  }
  clear();
}

search_result branch::search(boost::optional<search_entry_key> limit_key) {
  if (m_search_entries.empty()) {
    return search_result::STOP_SEARCHING;
  }

  if (!update_ploids_remaining()) {
    if (m_tracable_entry_count) {
      std::cout << "Clearing branch " << *this << " due to no ploids remaining\n";
    }
    clear();
    return search_result::STOP_SEARCHING;
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  while (!m_search_entries.empty()) {
    if (limit_key && m_search_entries.front()->get_key() < *limit_key) {
      break;
    }
    execute_one_search_internal();
  }
  auto end_time = std::chrono::high_resolution_clock::now();

  m_time_spent += end_time - start_time;
  return search_result::SEARCH_MORE;
}

void branch::execute_search_for_testing(branch_search_entry_ptr e) {
  update_ploids_remaining();
  execute_search_internal(std::move(e));
}

bool branch::trace_enabled_for_entry(const branch_search_entry* e) const {
  if (m_trace.empty()) {
    return false;
  }
  if (dynamic_cast<const rejoin_search_entry*>(e)) {
    return true;
  }
  if (trace_enabled(e->get_path())) {
    return true;
  }
  return false;
}

void branch::execute_search_internal(branch_search_entry_ptr e) {
  bool needs_trace = trace_enabled_for_entry(e.get());
  if (k_dbg || needs_trace) {
    std::cout << "Branch executing search (" << m_steps_left
              << " steps left) on: " << e->describe(this) << " (" << m_search_entries.size()
              << " searches left)\n";
  }
  unsigned orig_tracable_entry_count = m_tracable_entry_count;
  if (needs_trace) {
    CHECK_GT(m_tracable_entry_count, 0);
    --m_tracable_entry_count;
  }

  if (e->pair_match_count() < m_max_pair_match_count / 2) {
    if (needs_trace) {
      std::cout << "TRACE had too few pair matches\n";
    }
    e->notify_discard(this);
    return;
  }

  switch (e->search(this)) {
    case search_result::STOP_SEARCHING:
      if (needs_trace && m_tracable_entry_count < orig_tracable_entry_count) {
        std::cout << "LOST TRACE, returned STOP_SEARCHING\n";
      }
      break;
    case search_result::SEARCH_MORE:
      if (needs_trace && !trace_enabled_for_entry(e.get())) {
        std::cout << "DIVERGED FROM TRACE, wanting SEARCH_MORE: " << e->describe(this) << "\n";
      }
      add_search_entry(std::move(e));
      break;
    default:
      LOG(FATAL) << "Invalid search result";
  }
}

void branch::execute_one_search_for_testing() {
  update_ploids_remaining();
  execute_one_search_internal();
}

void branch::execute_one_search_internal() {
  CHECK(!m_search_entries.empty());

  std::pop_heap(m_search_entries.begin(), m_search_entries.end(), branch_search_entry_comparer());
  branch_search_entry_ptr e = std::move(m_search_entries.back());
  m_search_entries.pop_back();

  if (m_steps_left == 0) {
    if (trace_enabled_for_entry(e.get())) {
      std::cout << "DISCARD TRACED ENTRY: " << e->describe(this) << " due to out of steps\n";
      CHECK_GT(m_tracable_entry_count, 0);
      --m_tracable_entry_count;
    }
    e->notify_discard(this);
  } else {
    --m_steps_left;
    execute_search_internal(std::move(e));
  }
};

void branch::add_search_entry(branch_search_entry_ptr e) {
  if (trace_enabled_for_entry(e.get())) {
    ++m_tracable_entry_count;
    std::cout << "Added search entry that needs trace: " << e->describe(this) << "\n";
  }

  if (e->pair_match_count() > m_max_pair_match_count) {
    unsigned new_pair_matches = e->pair_match_count() - m_max_pair_match_count;
    m_steps_left += new_pair_matches * opts().bidir_branch_steps_per_pair;
    if (m_steps_left > opts().bidir_max_branch_steps) {
      m_steps_left = opts().bidir_max_branch_steps;
    }
    m_max_pair_match_count = e->pair_match_count();
  }

  m_search_entries.emplace_back(std::move(e));
  std::push_heap(m_search_entries.begin(), m_search_entries.end(), branch_search_entry_comparer());
}

void branch::check_path_invariants(const path& p) const {
  unsigned actual_anchor_len =
      push_view()->shared_ref_bases_to_left(right_push_view_offset() + p.anchor_len(), p.seq());
  CHECK_EQ(p.anchor_len(), actual_anchor_len)
      << "Right offset: " << right_push_view_offset() << " seq: " << p.seq()
      << "push search entry: " << *this << " scaffold: " << push_view()->opts().scaffold_name;

  if (push_view()->opts().bidir_validate_trace_state > 1) {
    p.check_invariants();
  }
}

void branch::clear() {
  // If we used all our steps getting here, free some up for potential
  // additional alleles.  But don't free up all of it...

  m_max_pair_match_count = 0;
  m_steps_left = opts().bidir_max_branch_steps;
  if (m_search_entries.empty()) {
    CHECK_EQ(0, m_tracable_entry_count);
  } else {
    std::vector<branch_search_entry_ptr> to_discard;
    std::swap(m_search_entries, to_discard);
    for (const auto& e : to_discard) {
      if (trace_enabled_for_entry(e.get())) {
        std::cout << "DISCARD TRACED ENTRY: " << e->describe(this) << " due to branch clearing\n";
        CHECK_GT(m_tracable_entry_count, 0);
        --m_tracable_entry_count;
      }

      e->notify_discard(this);
    }
    CHECK_EQ(m_tracable_entry_count, 0);
    CHECK(m_search_entries.empty());
  }
}

std::string branch_search_entry::describe(const branch* br) const {
  std::stringstream result;
  result << "BrSearch(" << m_key << "):" << describe_internal(br);
  return result.str();
}

bool branch::try_rejoin(aoffset_t outer_left_offset, dna_slice left_seq, const path& p,
                        unsigned pair_match_count) {
  if (k_dbg) {
    std::cout << "\nBranch considering rejoin at " << outer_left_offset << ":\n"
              << "Seq: " << left_seq << "\n"  //
              << "Path: " << p << "\n";
  }
  if (outer_left_offset < m_ref_remaining) {
    if (k_dbg) {
      std::cout << "Outer left " << outer_left_offset << " < ref remaining " << m_ref_remaining
                << "; cannot rejoin\n";
    }
    return false;
  }

  if (outer_left_offset >= right_push_view_offset()) {
    if (k_dbg) {
      std::cout << "Outer left " << outer_left_offset << " > right offset " << m_right_offset
                << "; cannot rejoin\n";
    }
    return false;
  }

  aoffset_t ref_distance = m_right_offset - outer_left_offset;
  if (ref_distance > aoffset_t(opts().read_ahead_distance)) {
    if (k_dbg) {
      std::cout << "ref distance " << ref_distance << " too far for readahead; cannot rejoin\n";
    }
    return false;
  }

  auto ext = push_view()->get_scaffold().split_extent_at(outer_left_offset);
  dna_slice ref_seq = ext.second;

  unsigned left_anchor_len = ref_seq.shared_prefix_length(left_seq);
  if (left_anchor_len == left_seq.size()) {
    left_anchor_len += ref_seq.subseq(left_anchor_len, ext.second.size() - left_anchor_len)
                           .shared_prefix_length(p.seq());
  }

  if (k_dbg) {
    std::cout << "Left anchor shares " << left_anchor_len << " bases with reference\nRef:\n"
              << ref_seq.subseq(0, left_anchor_len) << "\nLeft seq: " << left_seq << "\nPath: " << p
              << "\n";
    std::cout << "Ref continues: " << ref_seq.subseq(left_anchor_len, 30) << "\n";
  }

  if (left_anchor_len >= (left_seq.size() + p.size())) {
    if (k_dbg) {
      std::cout << "Left anchor is whole sequence; cannot rejoin\n";
    }
    return false;
  }

  CHECK_GE(left_anchor_len, opts().bidir_min_anchor_len)
      << "Left anchor only has " << left_anchor_len
      << " bases in common with reference; cannot rejoin\n";

  aoffset_t left_offset = outer_left_offset + left_anchor_len;

  path rejoin_path = p;
  rejoin_path.push_front_drop(left_seq);

  unsigned path_overlap = rejoin_path.path_overlap();
  path_overlap = std::min<unsigned>(path_overlap, left_anchor_len);

  if (k_dbg) {
    std::cout << "Rejoin path: " << rejoin_path << "\n";
  }

  if (left_offset >= right_push_view_offset() + p.anchor_len()) {
    return false;
  }

  std::unique_ptr<rejoin_search_entry> e = make_unique<rejoin_search_entry>(
      path_overlap, left_offset, left_anchor_len, std::move(rejoin_path), pair_match_count);
  if (opts().bidir_validate_trace_state) {
    e->check_invariants(this);
  }
  if (k_dbg) {
    std::cout << "Rejoin try successful; saving rejoin search entry:\n"
              << e->describe(this) << "\n";
  }

  add_search_entry(std::move(e));

  return true;
}

branch_search_entry::branch_search_entry(const search_entry_key& key) : m_key(key) {}

void branch::enable_trace(dna_slice seq) {
  CHECK_EQ(seq.rev_comp()[0].complement(), m_first_base);
  m_trace.emplace(seq.rev_comp());
}

bool branch::trace_enabled(dna_slice seq) const {
  if (m_trace.empty()) {
    return false;
  }

  dna_slice rc_seq = seq.rev_comp();

  auto it = m_trace.lower_bound(dna_sequence(rc_seq));
  if (it != m_trace.end()) {
    unsigned shared1 = rc_seq.shared_prefix_length(*it);
    CHECK_GE(shared1, 1);
    if (shared1 == rc_seq.size() || shared1 == it->size()) {
      return true;
    }
  }

  if (it == m_trace.begin()) {
    return false;
  }
  --it;

  unsigned shared2 = rc_seq.shared_prefix_length(*it);
  CHECK_GE(shared2, 1);
  if (shared2 == rc_seq.size() || shared2 == it->size()) {
    return true;
  }

  return false;
}

bool branch::trace_enabled(const path& p) const {
  if (m_trace.empty()) {
    return false;
  }

  dna_slice seq = p.seq();
  seq = seq.subseq(0, seq.size() - p.anchor_len());
  return trace_enabled(seq);
}

bool branch::explore(const seqset_range& r) { return m_explored.insert(r).second; }

std::ostream& operator<<(std::ostream& os, const branch& br) { return os << br.describe(); }

boost::optional<search_entry_key> branch::best_search_entry_key() const {
  if (m_search_entries.empty()) {
    return boost::none;
  }

  return m_search_entries.front()->get_key();
}

void branch::note_output(dna_slice seq) {
  if (opts().bidir_report_slow_branches) {
    m_outputs.emplace(seq);
  }
}

}  // namespace discovery
}  // namespace variants

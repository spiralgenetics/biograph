#include "modules/variants/discovery/state.h"
#include "modules/variants/discovery/path.h"
#include "modules/variants/discovery/rejoin.h"
#include "modules/variants/discovery/view.h"
#include "modules/variants/discovery/walk_ref.h"

#include "modules/bio_base/readmap.h"

static constexpr bool k_dbg = false;
namespace variants {
namespace discovery {

state::state(const assemble_options& options, pipeline_step_t owned_output)
    : m_options(options), m_owned_output(std::move(owned_output)) {
  std::tie(m_fwd_view, m_rev_view) = view_t::create_view_pair(*m_options.scaffold, this);
  CHECK(!m_fwd_view->is_rev_comp());
  CHECK(m_rev_view->is_rev_comp());
  m_ploids_remaining.insert(
      std::make_pair(interval_t(0, m_options.scaffold->end_pos()), m_options.bidir_max_ploids));
  m_output = m_owned_output.get();

  if (m_options.report_bidir_initialized_func) {
    m_options.report_bidir_initialized_func(this);
  }
}

state::~state() {}

bool state::explore(const seqset_range& r) { return m_explored.insert(r).second; }

std::vector<std::pair<search_entry_key, branch*>> state::all_search_entries() const {
  std::vector<std::pair<search_entry_key, branch*>> result;

  for (view_t* v : both_dirs()) {
    for (branch* br : v->branches()) {
      boost::optional<search_entry_key> best_key = br->best_search_entry_key();
      if (!best_key) {
        continue;
      }
      result.emplace_back(std::make_pair(*best_key, br));
    }
  }

  struct comparer {
    bool operator()(const std::pair<search_entry_key, branch*>& lhs,
                    const std::pair<search_entry_key, branch*>& rhs) const {
      return lhs.first < rhs.first;
    }
  };
  // Sort so best is at end.
  std::sort(result.begin(), result.end(), comparer());

  return result;
}

search_result state::execute_one_search() {
  if (m_options.bidir_validate_trace_state > 1) {
    check_invariants();
  }

  std::vector<std::pair<search_entry_key, branch*>> search_entries = all_search_entries();

  if (search_entries.empty()) {
    return search_result::STOP_SEARCHING;
  }

  while (!search_entries.empty()) {
    branch* br = search_entries.back().second;
    search_entries.pop_back();

    if (search_entries.empty()) {
      br->search(boost::none /* search until done */);
      break;
    }

    search_entry_key next_best = search_entries.back().first;
    br->search(next_best);
  }

  return search_result::SEARCH_MORE;
}

void state::output_assembly(assembly_ptr a, bool uses_ploid) {
  if (m_options.bidir_validate_trace_state) {
    check_assembly(*a, "discovery output assembly");
  }

  if (uses_ploid) {
    m_ploids_remaining -=
        std::make_pair(interval_t(a->left_offset + aoffset_t((a->left_anchor_len + 1) / 2),
                                  a->right_offset - aoffset_t((a->right_anchor_len + 1) / 2)),
                       ploids_remaining_counter(1));
  }

  m_output->add(std::move(a));
}

void state::add_reference(aoffset_t start, aoffset_t limit) {
  for (view_t* v : both_dirs()) {
    walk_ref_t wr(v);

    aoffset_t dir_start, dir_limit;

    if (v->is_rev_comp()) {
      dir_start = v->reverse_offset(limit);
      dir_limit = v->reverse_offset(start);
    } else {
      dir_start = start;
      dir_limit = limit;
    }
    wr.walk_ref(dir_start, dir_limit);
    if (opts().bidir_validate_trace_state) {
      wr.check_invariants();
    }
    wr.init_pairs_and_push();
  }
}

void state::assemble(assemble_pipeline_interface* output, progress_handler_t progress) {
  CHECK(!m_output || !output);

  if (output) {
    m_output = output;
  }
  while (execute_one_search() == search_result::SEARCH_MORE) {
  }
  progress(1);
  if (output) {
    m_output = nullptr;
  }

  if (opts().bidir_report_slow_branches) {
    show_longest_branches();
  }
}

void state::show_longest_branches() {
  std::vector<branch*> branches;
  for (view_t* v : both_dirs()) {
    for (branch* br : v->branches()) {
      branches.emplace_back(br);
    }
  }

  if (branches.empty()) {
    return;
  }

  struct comparer {
    bool operator()(branch* a, branch* b) const { return a->time_spent() < b->time_spent(); }
  };

  auto max_it = std::max_element(branches.begin(), branches.end(), comparer());
  CHECK(max_it != branches.end());

  static std::mutex g_mu;
  std::lock_guard<std::mutex> l(g_mu);
  branch* br = *max_it;
  static std::chrono::nanoseconds max_time_spent{0};

  if (br->time_spent() < max_time_spent) {
    return;
  }

  const auto& outs = br->outputs();
  if (!outs.empty()) {
    return;
  }
  max_time_spent = br->time_spent();
  std::cout << "New longest branch " << *br << " took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(max_time_spent).count()
            << " ms and produced " << outs.size() << " outputs:\n";
  for (const auto& out : outs) {
    std::cout << "  " << out << "\n";
  }
}

offset_info state::get_offset_info(aoffset_t offset, bool fwd) const {
  if (k_dbg) {
    std::cout << "Attempting to get offset info at " << offset << " fwd=" << fwd << "\n";
    std::cout << "Table is: " << m_ploids_remaining << "\n";
  }

  offset_info result;
  result.ploids_remaining = 0;

  auto it = m_ploids_remaining.find(offset);
  if (it == m_ploids_remaining.end()) {
    if (k_dbg) {
      std::cout << "Could not find offset in table\n";
    }
    result.ploids_remaining = 0;
    return result;
  }
  if (k_dbg) {
    std::cout << "Found interval " << it->first << " With ploids remaining=" << it->second << "\n";
  }
  CHECK_GT(it->second, 0) << it->first;
  result.ploids_remaining = it->second;

  if (fwd) {
    result.ref_remaining_limit = it->first.lower() - 1;
    if (k_dbg) {
      std::cout << "Starting forward trace at " << it->first << ": " << result << "\n";
    }
    while (it != m_ploids_remaining.end() && it->first.lower() - 1 == result.ref_remaining_limit) {
      result.ref_remaining_limit = it->first.upper();
      if (k_dbg) {
        std::cout << "Continuing forward trace at " << it->first << ": " << result << "\n";
      }
      CHECK_GT(it->second, 0) << it->first;
      ++it;
    }
    if (k_dbg) {
      std::cout << "Done forward trace\n";
    }
  } else {
    result.ref_remaining_limit = it->first.upper() + 1;
    if (k_dbg) {
      std::cout << "Starting reverse trace at " << it->first << ": " << result << "\n";
    }
    for (;;) {
      if (k_dbg) {
        std::cout << "Continuing reverse trace at " << it->first << ": " << result << "\n";
      }
      if (result.ref_remaining_limit == it->first.upper() + 1) {
        CHECK_GT(it->second, 0) << it->first;
        result.ref_remaining_limit = it->first.lower();
      } else {
        break;
      }
      if (it != m_ploids_remaining.begin()) {
        --it;
      } else {
        break;
      }
    }
  }

  if (k_dbg) {
    std::cout << "Result offset info: " << result << "\n";
  }
  return result;
}

void state::check_invariants() const {
  m_fwd_view->check_invariants();
  m_rev_view->check_invariants();

  for (const auto& elem : m_ploids_remaining) {
    CHECK_LE(elem.second, m_options.bidir_max_ploids);
    CHECK_GT(elem.second, 0) << elem.first;
  }
}

std::ostream& operator<<(std::ostream& os, const offset_info& oi) {
  return os << "OffsetInfo(ploids_remaining=" << oi.ploids_remaining
            << ",ref_remaining_limit=" << oi.ref_remaining_limit << ")";
}

void state::discard_search_entries() {
  for (view_t* v : both_dirs()) {
    v->discard_search_entries();
  }
}

void state::add_trace_for_variant(aoffset_t left_offset, aoffset_t right_offset, dna_slice seq) {
  auto left_ext = fwd_view()->get_scaffold().split_extent_at(left_offset);
  auto right_ext = fwd_view()->get_scaffold().split_extent_at(right_offset);

  // Calculate full sequence including variant extending to the right starting at left_offset.
  dna_sequence left;
  left += seq;
  left += right_ext.second;

  // Calculate full sequence including variant extending to the left starting at right_offset.
  dna_sequence rc_right;
  rc_right += seq.rev_comp();
  rc_right += left_ext.first.rev_comp();

  unsigned left_shared = left.shared_prefix_length(left_ext.second);
  unsigned right_shared = rc_right.shared_prefix_length(right_ext.first.rev_comp());

  CHECK_LT(left_shared, left.size());
  CHECK_LE(left_shared, left_ext.second.size());

  CHECK_LT(right_shared, rc_right.size());
  CHECK_LE(right_shared, right_ext.first.rev_comp().size());

  aoffset_t adjusted_left = left_offset + left_shared;
  dna_base left_branch_base = left[left_shared];

  aoffset_t adjusted_right = right_offset - right_shared;
  dna_base rc_right_branch_base = rc_right[right_shared];

  fwd_view()
      ->get_branch(rc_right_branch_base.complement(), adjusted_right)
      ->enable_trace(
          dna_slice(rc_right).subseq(right_shared, rc_right.size() - right_shared).rev_comp());
  rev_view()
      ->get_branch(left_branch_base.complement(), rev_view()->reverse_offset(adjusted_left))
      ->enable_trace(dna_slice(left).subseq(left_shared, left.size() - left_shared).rev_comp());
}

std::ostream& operator<<(std::ostream& os, const search_entry_key& e) {
  os << "ol=" << e.path_overlap << ",";
  switch (e.priority) {
    case search_priority::PUSH:
      os << "PUSH";
      break;
    case search_priority::POP:
      os << "POP";
      break;
    case search_priority::REJOIN:
      os << "REJOIN";
      break;
  }
  os << ",p=" << e.pair_match_count;
  return os;
}

bool search_entry_key::operator<(const search_entry_key& bkey) const {
  const search_entry_key& akey = *this;

  // Search more pairs matched first.
  if (akey.pair_match_count != bkey.pair_match_count) {
    return bkey.pair_match_count > akey.pair_match_count;
  }

  // Try things with the best minimum overlap first.
  if (akey.path_overlap != bkey.path_overlap) {
    return bkey.path_overlap > akey.path_overlap;
  }

  // Search better average overlap first.
  unsigned a_avg_overlap = akey.tot_overlap * bkey.tot_overlap_bases;
  unsigned b_avg_overlap = bkey.tot_overlap * akey.tot_overlap_bases;
  if (a_avg_overlap != b_avg_overlap) {
    return b_avg_overlap > a_avg_overlap;
  }

  // Try things with a better priority first
  if (akey.priority != bkey.priority) {
    return bkey.priority > akey.priority;
  }

  return false;
};

search_entry_key::search_entry_key(search_priority new_priority, const path& p,
                                   unsigned new_pair_match_count)
    : path_overlap(p.path_overlap()),
      tot_overlap(p.tot_overlap()),
      tot_overlap_bases(p.tot_overlap_bases()),
      priority(new_priority),
      pair_match_count(new_pair_match_count) {}

}  // namespace discovery
}  // namespace variants

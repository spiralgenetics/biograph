#pragma once

#include "modules/bio_base/readmap.h"
#include "modules/variants/discovery/seqset_range_table.h"
#include "modules/variants/discovery/state.h"
#include "modules/variants/discovery/view.h"

namespace variants {
namespace discovery {

class branch_search_entry {
 public:
  virtual ~branch_search_entry() = default;

  search_result search(branch* br);

  virtual void check_invariants(const branch* br) const = 0;

  std::string describe(const branch* br) const;

  unsigned path_overlap() const { return get_key().path_overlap; }
  search_priority priority() const { return get_key().priority; }
  unsigned pair_match_count() const { return get_key().pair_match_count; }

  virtual const path& get_path() const = 0;

  virtual void notify_discard(branch* br) {}

  const search_entry_key& get_key() const { return m_key; }

 protected:
  branch_search_entry() = delete;
  branch_search_entry(const search_entry_key& new_key);
  branch_search_entry& operator=(const branch_search_entry&) = delete;

  virtual search_result search_internal(branch* br) = 0;
  virtual std::string describe_internal(const branch* br) const = 0;
  virtual unsigned cur_overlap() const { return 0; }

  search_entry_key m_key;
};

using branch_search_entry_ptr = std::unique_ptr<branch_search_entry>;

// Represents a branch off of reference.  Keeps a queue of search entries to process.
class branch {
 public:
  // Creates a new branch extending to the left starting at right_offset.  The first base in the
  // branch is first_base.
  branch(view_t* push_view, dna_base first_base, aoffset_t right_offset);

  // Right offset of the branch relative to push_view
  aoffset_t right_push_view_offset() const { return m_right_offset; }
  // Left offset of the branch relative to pop_view.
  aoffset_t left_pop_view_offset() const { return m_push_view->reverse_offset(m_right_offset); }

  aoffset_t pop_view_farthest_right_offset() const {
    if (m_ref_remaining == std::numeric_limits<aoffset_t>::max()) {
      return std::numeric_limits<aoffset_t>::min();
    } else {
      return m_push_view->reverse_offset(m_ref_remaining);
    }
  }
  aoffset_t push_view_farthest_left_offset() const { return m_ref_remaining; }

  dna_base first_base() const { return m_first_base; }
  bool empty() const { return m_search_entries.empty(); }

  // Search this branch for a while.  If limit_key isn't boost::none,
  // stop searching if we become worse than limit_key.
  search_result search(boost::optional<search_entry_key> limit_key);

  view_t* push_view() const;
  view_t* pop_view() const;
  state* get_state() const;

  void add_search_entry(branch_search_entry_ptr e);

  void check_invariants() const;
  std::string describe() const;

  boost::optional<search_entry_key> best_search_entry_key() const;

  void check_path_invariants(const path& p) const;

  void clear();

  void notify_rejoin(branch* other_br, branch_search_entry* other_e);
  const assemble_options& opts() const { return m_push_view->opts(); }

  // Adds a rejoin to the search queue.  Returns true if successful.
  // outer_left_offset is at the beginning of left_seq, and is to the
  // left of the reference anchor on the left.
  //
  // The rejoin assembly is left_seq + p.seq().
  bool try_rejoin(aoffset_t outer_left_offset, dna_slice left_seq, const path& p,
                  unsigned pair_match_count);

  // Debugging access:
  void enable_trace(dna_slice seq);
  bool trace_enabled(dna_slice seq) const;
  bool trace_enabled(const path& p) const;
  bool trace_enabled_for_entry(const branch_search_entry* e) const;
  bool any_trace_enabled() const { return !m_trace.empty(); }

  void note_output(dna_slice seq);
  const std::set<dna_sequence>& outputs() const { return m_outputs; }

  // Returns false if r has been previously explored.
  bool explore(const seqset_range& r);

  const std::chrono::nanoseconds& time_spent() const { return m_time_spent; }

 private:
  friend class discovery_test;
  void execute_one_search_for_testing();
  void execute_search_for_testing(branch_search_entry_ptr);

  void execute_one_search_internal();
  void execute_search_internal(branch_search_entry_ptr);
  // Returns false if no ploids are available.
  bool update_ploids_remaining();

  view_t* const m_push_view;
  dna_base const m_first_base;
  aoffset_t const m_right_offset;

  aoffset_t m_ref_remaining = std::numeric_limits<aoffset_t>::max();

  std::vector<branch_search_entry_ptr> m_search_entries;

  // Seqset entries that have already been explored, to avoid duplicating work.
  std::unordered_set<seqset_range, seqset_range_hash> m_explored;

  // Sequences we should output debugging traces for
  std::set<dna_sequence> m_trace;

  // Number of entries in the search queue that are tracable.
  unsigned m_tracable_entry_count = 0;

  unsigned m_steps_left = 0;

  // Maximum number of pair matches we've seen on any path off this branch.
  unsigned m_max_pair_match_count = 0;

  std::set<dna_sequence> m_outputs;

  // Total time spent searching through this branch.
  std::chrono::nanoseconds m_time_spent{0};
};

std::ostream& operator<<(std::ostream& os, const branch& br);

}  // namespace discovery
}  // namespace variants

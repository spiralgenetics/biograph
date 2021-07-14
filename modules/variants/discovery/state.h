#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval_set.hpp>

#include <map>
#include <memory>
#include <vector>

namespace variants {
namespace discovery {

class ploids_remaining_counter {
 public:
  ploids_remaining_counter() = default;
  ploids_remaining_counter(int val) : m_val(val) {}

  explicit operator bool() const { return m_val > 0; }
  operator int() const { return m_val; }
  ploids_remaining_counter& operator-=(const ploids_remaining_counter& to_add) {
    m_val -= to_add.m_val;
    if (m_val < 0) {
      m_val = 0;
    }
    return *this;
  }

 private:
  int m_val = 0;
};

using interval_t = boost::icl::closed_interval<aoffset_t>;
using interval_set_t = boost::icl::interval_set<aoffset_t, std::less, interval_t>;
using ploids_remaining_t =
    boost::icl::interval_map<aoffset_t, ploids_remaining_counter, boost::icl::partial_absorber,
                             std::less, boost::icl::inplace_plus, boost::icl::inter_section,
                             interval_t>;

enum class search_result { STOP_SEARCHING, SEARCH_MORE };
enum class search_priority { POP, PUSH, REJOIN };
class path;
struct search_entry_key {
  search_entry_key() = delete;
  search_entry_key(search_priority new_priority, const path& p, unsigned new_pair_match_count);

  // Minimum overlap along the path
  unsigned path_overlap;

  // Total overlap and base count so we can calculate the average.
  unsigned tot_overlap;
  unsigned tot_overlap_bases;

  search_priority priority;
  unsigned pair_match_count;

  // Orders search entry keys from worst to best.
  bool operator<(const search_entry_key& rhs) const;
  bool operator>(const search_entry_key& rhs) const { return rhs < *this; }
};

std::ostream& operator<<(std::ostream& os, const search_entry_key& key);
void PrintTo(const interval_t& intv, std::ostream* os);

struct offset_info {
  // Number of ploids at this reference offset that haven't been output.
  int ploids_remaining = std::numeric_limits<int>::min();

  // First reference position (either going forward or reverse,
  // depending on the arugment to get_offset_info) that has no ploids
  // available for output.
  aoffset_t ref_remaining_limit = std::numeric_limits<aoffset_t>::max();
};

std::ostream& operator<<(std::ostream& os, const offset_info& oi);

// Current state of discovery.
class view_t;
class branch;
class state {
 public:
  state(const assemble_options& options, pipeline_step_t output = nullptr);
  state() = delete;
  state(const state&) = delete;
  ~state();

  void add_reference(aoffset_t ref_start, aoffset_t ref_limit);
  void assemble(assemble_pipeline_interface* output = nullptr,
                progress_handler_t progress = null_progress_handler);

  search_result execute_one_search();

  // If uses_ploid is true, counts the emitted assembly against the ploids avaliable.
  // bidir_tracer_emit_all_rejoins may call this with uses_ploid=false.
  void output_assembly(assembly_ptr a, bool uses_ploid);

  // Returns false if r has been previously explored.
  bool explore(const seqset_range& r);

  // If "fwd" is true, ref_remaining_limit will be >= position.
  // If "fwd" is false, ref_remaining_limit will be <= position.
  // If ploids_remaining is <= 0, ref_remaining_limit is undefined.
  offset_info get_offset_info(aoffset_t position, bool fwd) const;

  std::array<view_t*, 2> both_dirs() const { return {{m_fwd_view.get(), m_rev_view.get()}}; }

  const assemble_options& opts() const { return m_options; }

  // Testing access:
  void discard_search_entries();
  view_t* fwd_view() const { return m_fwd_view.get(); }
  view_t* rev_view() const { return m_rev_view.get(); }

  // Make sure all invariants are fulfilled.
  void check_invariants() const;

  // Debugging access:
  void add_trace_for_variant(aoffset_t left_offset, aoffset_t right_offset, dna_slice seq);

 private:
  friend class discovery_test;
  friend class state_test;

  void show_longest_branches();
  std::vector<std::pair<search_entry_key, branch*>> all_search_entries() const;

  // Number of ploids remaining to output.  Initially this gets set
  // to m_options.bidir_max_ploids, and gets decreased as things get output.
  ploids_remaining_t m_ploids_remaining;

  // Seqset entries that have already been explored, to avoid duplicating work.
  std::unordered_set<seqset_range, seqset_range_hash> m_explored;

  std::unique_ptr<view_t> m_fwd_view;
  std::unique_ptr<view_t> m_rev_view;

  assemble_options m_options;
  assemble_pipeline_interface* m_output = nullptr;
  pipeline_step_t m_owned_output;
};

}  // namespace discovery
}  // namespace variants

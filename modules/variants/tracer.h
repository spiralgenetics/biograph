#pragma once

#include <boost/optional.hpp>
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/variants/assemble.h"
#include "modules/variants/assemble.h"
#include "modules/variants/ref_map.h"
#include "modules/variants/scaffold.h"

namespace variants {

class tracer {
 public:
  static const char k_tracer_name[];

  tracer(const assemble_options& options);

  tracer() = delete;
  tracer(const tracer&) = delete;
  ~tracer();

  assemble_stats assemble(assemble_pipeline_interface* output,
                          progress_handler_t progress = null_progress_handler);
  assemble_stats assemble(aoffset_t start_offset, aoffset_t limit_offset,
                          assemble_pipeline_interface* output,
                          progress_handler_t = null_progress_handler);

  void output_path_debug_dot(std::ostream& os) const;

 private:
  class path_storage;
  struct path_debug_info;
  struct path;
  struct allele_entry;
  struct next_path;
  struct rejoin;
  struct dead_end_sorter;

  struct next_path_comparer;
  struct rejoin_comparer {
    // Compares rejoin candidates.  Returns true if rhs is better than p.
    bool operator()(const rejoin& p, const rejoin& rhs) const;
  };

  void expand_next_path(next_path np);
  void add_base_to_next_path(next_path& np, dna_base b);
  void add_dead_end_rejoin(const next_path& np);
  bool add_rejoins(const next_path& np);
  void add_pairs_to_next_path(next_path& np);

  void output_rejoins(assemble_pipeline_interface* output);
  void output_dead_ends(unsigned max_to_output, assemble_pipeline_interface* output);

  void output_assembly(const rejoin& r,
                       assemble_pipeline_interface* output) const;

  void push_next_path(next_path np);
  next_path pop_next_path();
  bool path_has_read_in_range(const path* start_path, uint32_t read_id,
                              aoffset_t start, aoffset_t limit) const;
  void push_rejoin(rejoin r);
  bool is_read(const seqset_range& r) const;
  std::string next_path_to_string(const next_path& np) const;
  std::string rejoin_to_string(const rejoin& r) const;
  std::string path_to_string(const path* p) const;

  void trace();
  void trace_next_paths();
  void skip_ahead_to(aoffset_t offset);
  void advance_read_ahead_to(aoffset_t offset);
  void advance_trail_behind_to(aoffset_t offset);

  bool has_seqset_id_in_range(uint64_t seqset_id, aoffset_t start_offset,
                              aoffset_t limit_offset) const;
  bool has_rc_mate_in_range(uint32_t read_id, aoffset_t start_offset,
                            aoffset_t limit_offset) const;

  void advance_position_entry_index();

  // Configuration:
  const std::vector<scaffold::extent>* m_ref_parts = nullptr;
  const assemble_options m_options;
  const seqset* const m_seqset;
  const readmap* const m_readmap;
  const reference* const m_ref;
  const ref_map* const m_rmap;

  // State:
  aoffset_t m_start_offset = 0;
  aoffset_t m_limit_offset = std::numeric_limits<aoffset_t>::max();
  size_t m_position_entry_index = 0;

  // Reads:
  std::deque<std::pair<aoffset_t /* offset */, uint64_t /* seqset id */>>
      m_position_entries;
  std::unordered_multimap<uint64_t /* seqset_id */, aoffset_t /* offset */, unsalted_hash>
      m_entry_positions;
  std::unordered_multimap<uint32_t /* read_id */, aoffset_t /* offset */, unsalted_hash>
      m_rc_mate_read_positions;

  // Readahead:
  scaffold::iterator m_read_ahead_it;
  seqset_range m_read_ahead_range;
  aoffset_t m_read_ahead_offset = std::numeric_limits<aoffset_t>::min();

  // Current trace state:
  std::unique_ptr<path_storage> m_trace_path_storage;
  path* m_prev_path = nullptr;
  aoffset_t m_cur_offset = 0;
  seqset_range m_cur_range;
  dna_sequence m_cur_left_anchor;
  size_t m_cur_ref_ambiguous_bases = 0;
  std::vector<next_path> m_next_paths;
  std::multiset<rejoin, rejoin_comparer> m_rejoin_paths;
  std::set<const path*> m_dead_end_rejoins;
  // Number of outputs so far on this trace.
  size_t m_trace_outputs = 0;
  // Search steps so far
  size_t m_search_step_count = 0;
  size_t m_ambiguous_search_step_count = 0;

  mutable std::map<const path*, path_debug_info> m_path_debugs;
  mutable std::unique_ptr<path_storage> m_debug_path_storage;

  mutable assemble_stats m_stats;
};

struct tracer::path {
  const path* prev = nullptr;

  acost_t cost = 0;
  uint8_t min_overlap = 0;
  seqset_range range;
  dna_sequence seq;
  mutable bool part_of_assembly = false;

  struct seen_pair {
    uint32_t read_id = std::numeric_limits<uint32_t>::max();
    aoffset_t offset = std::numeric_limits<uint32_t>::max();

    bool operator==(const seen_pair& rhs) const {
      return read_id == rhs.read_id && offset == rhs.offset;
    }
    bool operator<(const seen_pair& rhs) const {
      if (read_id != rhs.read_id) {
        return read_id < rhs.read_id;
      }
      return offset < rhs.offset;
    }
  };
  std::vector<seen_pair> seen_pairs;
  mutable std::vector<uint32_t> seen_read_ids;
};

struct tracer::next_path {
  // Guaranteed that nothing else has a pointer to this particular
  // path until the next_path entry has been expanded.
  path* new_path = nullptr;

  unsigned path_bases = 0;

  unsigned pushed_since_read = 0;
  unsigned pushed_since_pair = 0;
  unsigned ambiguous_bases = 0;
  unsigned branch_count_since_pair = 0;
  unsigned max_between_pairs = 0;
  unsigned pairs_used = 0;
  unsigned num_reads = 0;

  uint64_t loop_check_seqset_id = std::numeric_limits<uint64_t>::max();
  unsigned reads_until_loop_check = 0;
};

struct tracer::rejoin {
  const path* p = nullptr;

  acost_t rejoin_cost = 0;
  aoffset_t right_offset = std::numeric_limits<aoffset_t>::max();
  int right_anchor_len = 0;
  bool anchor_drop = false;
};

}  // namespace variants

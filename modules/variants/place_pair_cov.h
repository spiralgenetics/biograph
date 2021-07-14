#pragma once

// Places pairs based on closeness to the specified ideal distance.
//
// Any reads in pair_coverage that share the same read_id will be placed once, based on
// closeness to the ideal distance.
//
// Incoming assemblies are assumed to already have valid pairs calculated.
#include "modules/variants/assemble.h"
#include "modules/variants/dist_set.h"

namespace variants {

struct place_pair_options {
  // Ideal distance between ends of pairs
  aoffset_t ideal_pair_distance = std::numeric_limits<aoffset_t>::max();

  // Maximum number of equally good places to place a pair.  Any pairs
  // which have more than this number equally good places will be
  // discarded.
  unsigned max_ambig = 15;
};

class place_pair_cov : public sorted_output_pipeline_step {
 public:
  using placer_inspector_t = std::function<void(place_pair_cov*)>;
  place_pair_cov(const assemble_options& opts, const place_pair_options& popts,
                 pipeline_step_t output);
  ~place_pair_cov();

  void on_assembly(assembly_ptr a) override;
  void flush() override;

  // For testing access:
  absl::btree_set<aoffset_t> testing_distances_between(aoffset_t left_offset,
                                                       aoffset_t right_offset) const;
  void dump_state(const std::string& where);

  void testing_set_inspector(placer_inspector_t inspect_for_testing) {
    m_inspect_for_testing = inspect_for_testing;
  }

 private:
  using dist_set_t = absl::btree_set<aoffset_t>;
  using dists_t = absl::btree_map<aoffset_t /* left offset */, dist_set>;
  std::string dump_dist_table(const dists_t& dt) const;
  struct dist_set_formatter;
  struct trace_state {
    trace_state(size_t trace_id_arg, assembly* a_arg) : trace_id(trace_id_arg), a(a_arg) {}

    size_t const trace_id;
    assembly* const a;

    read_coverage_set filtered_coverage;
    std::vector<trace_state*> left_edges;
    std::vector<trace_state*> right_edges;

    // Round robin indexes for spreading out alignments.
    size_t left_edge_rr_idx = 0;
    size_t right_edge_rr_idx = 0;
  };
  friend std::ostream& operator<<(std::ostream& os, const trace_state& st);

  struct align_part {
    trace_state* st = nullptr;
    int offset = std::numeric_limits<aoffset_t>::max();

    bool operator==(const align_part& rhs) const { return st == rhs.st && offset == rhs.offset; }
  };
  friend std::ostream& operator<<(std::ostream& os, const align_part& part);

  struct anchor {
    // Assembly which this anchor occurs in.
    trace_state* st = nullptr;

    // Offset from beginning of assembly of the anchored read
    aoffset_t offset = std::numeric_limits<aoffset_t>::max();

    // Ordering that's deterministic and doesn't depend on where the
    // "trace_state" pointers are allocated.
    bool operator<(const anchor& rhs) const {
      const auto& lhs = *this;
      if (lhs.st != rhs.st) {
        return lhs.st->trace_id < rhs.st->trace_id;
      }
      return lhs.offset < rhs.offset;
    }
  };

  struct align {
    std::vector<align_part> parts;

    // Reference anchor
    aoffset_t ref_anchor = std::numeric_limits<aoffset_t>::max();
    trace_state* ref_anchor_st = nullptr;
  };

  struct pair_align_metric {
    int dist_from_ideal = std::numeric_limits<aoffset_t>::max();
  };
  friend std::ostream& operator<<(std::ostream& os, const pair_align_metric& metric);
  struct gather_anchors {
    pair_align_metric best_metric;
    std::vector<std::pair<const anchor*, const anchor*>> best_pairs;
  };

  using anchors_t = absl::btree_set<anchor>;
  struct read_info {
    // Right anchors of left mates
    anchors_t read_states;

    // Left anchors of right mates
    anchors_t rc_mate_states;
  };
  friend std::ostream& operator<<(std::ostream& os, const read_info& st);

  void place();

  void calc_dists();
  void save_reads(assembly* a, trace_state* st);
  void init_edges();
  void filter_reads();
  void save_filtered_reads();

  trace_state* assembly_to_state(assembly* a);

  void add_dists(const dists_t& old, aoffset_t distance, dists_t& result);

  void place_and_filter(uint32_t read_id, const read_info& ri);
  bool state_has_read(trace_state* st, uint32_t read_id, aoffset_t offset, int read_len) const;

  // Returns true if any reads were successfully propagated.
  bool propagate_read(uint32_t read_id, trace_state* st, aoffset_t offset, aoffset_t read_len,
                      bool prop_left /* true if propagating left, false if propagating right */,
                      std::vector<align>& align_out, align& so_far);

  void gather_anchor(gather_anchors& anchors, const pair_align_metric& metric, const anchor* left,
                     const anchor* right, bool brief_dbg);
  void save_align(const align& aln, uint32_t read_id, int read_len);
  void dump_align(const align& aln, uint32_t read_id, int read_len) const;
  bool metric_better(const pair_align_metric& lhs, const pair_align_metric& rhs) const;

  void debug_assembly(const assembly& a);
  void maybe_debug_assembly(const assembly& a);

  size_t m_next_trace_id = 0;
  std::vector<assembly_ptr> m_block;

  absl::flat_hash_map<assembly*, std::unique_ptr<trace_state>> m_asm_to_state;

  absl::btree_map<uint32_t /* left mate read id */, read_info> m_reads;

  // Distances to each reference location from previous locations.
  // The distances from reference location 10 to 100 are in m_dists[100][10].
  using dists_table_t = absl::btree_map<aoffset_t /* right offset */, dists_t>;
  dists_table_t m_dists;

  // Maximum distance we need to keep track of between any two reference points.  Any distances
  // loner than this are irrelevant, e.g. m_max_dist + 2 * read_len <= max_pair_distance
  aoffset_t m_max_dist = 0;

  // Maximum ideal distance we need to keep track of between any two
  // reference points.  We only need to track a single distance longer than this.
  aoffset_t m_max_ideal_dist = 0;

  // Index to use for round robining between alignments.
  size_t m_rr_idx = 0;

  assemble_options m_opts;
  place_pair_options m_popts;

  placer_inspector_t m_inspect_for_testing;

  // Temporaries to use to avoid per-read allocations
  std::vector<anchor> m_lefts, m_rights;
  align m_so_far;

  static absl::btree_set<uint32_t> g_debug_read_ids;
  static absl::btree_set<size_t> g_debug_assembly_ids;
};

}  // namespace variants

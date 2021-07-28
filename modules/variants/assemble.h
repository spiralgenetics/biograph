#pragma once

#include <boost/any.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/optional.hpp>
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/io/autostats.h"
#include "modules/io/ref_count.h"
#include "modules/variants/read_set.h"
#include "modules/variants/ref_map.h"

namespace variants {
namespace discovery {
class state;
}

// Assembly offset.
using aoffset_t = int32_t;

// Assembly search cost.
using acost_t = int64_t;

struct assembly;
// TODO(nils): Make sure we don't lose asemblies by disallowing explicit copy and delete here.
using assembly_ptr =
    explicit_shared_ptr<assembly, true /* atomic */, true /* allow implicit copy */,
                        true /* allow implicit delete */>;

// Version of boost::optional that autoconverts into aoffset or bool,
// and throws an exception if the offset is not present.
//
// Just a wrapper around boost::optional so we don't have to do stuff
// like "*a->left_offset" everywhere.
class optional_aoffset {
 public:
  optional_aoffset() : m_val(boost::none) {}
  optional_aoffset(const aoffset_t& val) : m_val(val) {}
  optional_aoffset(const optional_aoffset&) = default;

  optional_aoffset& operator=(const optional_aoffset& val) {
    m_val = val.m_val;
    return *this;
  }
  optional_aoffset& operator=(const aoffset_t& val) {
    m_val = val;
    return *this;
  }

  optional_aoffset& operator+=(aoffset_t rhs) {
    mut_val() += rhs;
    return *this;
  }
  optional_aoffset& operator-=(aoffset_t rhs) {
    mut_val() -= rhs;
    return *this;
  }
  aoffset_t operator--(int) { return mut_val()--; }
  aoffset_t operator++(int) { return mut_val()--; }
  optional_aoffset& operator--() {
    --mut_val();
    return *this;
  }
  optional_aoffset& operator++() {
    ++mut_val();
    return *this;
  }

  aoffset_t operator+(aoffset_t rhs) const { return get_val() + rhs; }
  aoffset_t operator+(size_t rhs) const { return get_val() + rhs; }
  aoffset_t operator-(aoffset_t rhs) const { return get_val() - rhs; }
  aoffset_t operator-(size_t rhs) const { return get_val() - aoffset_t(rhs); }

  bool operator==(const optional_aoffset& rhs) const { return m_val == rhs.m_val; }
  bool operator!=(const optional_aoffset& rhs) const { return m_val != rhs.m_val; }
  bool operator==(aoffset_t rhs) const { return get_val() == rhs; }
  bool operator!=(aoffset_t rhs) const { return get_val() != rhs; }
  bool operator<(aoffset_t rhs) const { return get_val() < rhs; }
  bool operator<=(aoffset_t rhs) const { return get_val() <= rhs; }
  bool operator>(aoffset_t rhs) const { return get_val() > rhs; }
  bool operator>=(aoffset_t rhs) const { return get_val() >= rhs; }
  bool operator==(size_t rhs) const { return get_val() == aoffset_t(rhs); }

  explicit operator bool() const { return m_val ? true : false; }
  operator const aoffset_t&() const { return get_val(); }
  friend std::ostream& operator<<(std::ostream& os, const optional_aoffset& ao) {
    if (ao) {
      return os << *ao.m_val;
    } else {
      return os << "(unanchored)";
    }
  }

  static const optional_aoffset none;

 private:
  const aoffset_t& get_val() const {
    // In debugging mode, make sure we get a stack trace:
    DCHECK(m_val);

    if (!m_val) {
      throw(std::runtime_error("Missing assembly offset"));
    }
    return *m_val;
  };
  aoffset_t& mut_val() {
    // In debugging mode, make sure we get a stack trace:
    DCHECK(m_val);

    if (!m_val) {
      throw(std::runtime_error("Missing assembly offset"));
    }
    return *m_val;
  };

  boost::optional<aoffset_t> m_val;
};

optional_aoffset min(const optional_aoffset& lhs, const optional_aoffset& rhs);
optional_aoffset max(const optional_aoffset& lhs, const optional_aoffset& rhs);

struct aligned_var {
  aoffset_t left_offset;
  aoffset_t right_offset;
  dna_sequence seq;

  // Result from genotyping:
  int max_alt_depth = 0;

  bool operator<(const aligned_var& rhs) const {
    if (left_offset != rhs.left_offset) {
      return left_offset < rhs.left_offset;
    }
    if (right_offset != rhs.right_offset) {
      return right_offset < rhs.right_offset;
    }
    if (seq != rhs.seq) {
      return seq < rhs.seq;
    }
    return false;
  }
  bool operator==(const aligned_var& rhs) const {
    return left_offset == rhs.left_offset && right_offset == rhs.right_offset && seq == rhs.seq;
  }
  bool operator!=(const aligned_var& rhs) const { return !(*this == rhs); }

  bool empty() const { return left_offset == right_offset && seq.size() == 0; }

  friend std::ostream& operator<<(std::ostream& os, const aligned_var& v);
};

struct align_count_t {
  // Sum of lengths of distinct reads.
  size_t local_read_lens = 0;

  // Sum of lengths of first alignment of reads in this assembly.
  size_t local_aligned_bases = 0;

  // Total aligned bases across all alignments for all reads aligned
  // in this assembly.
  size_t tot_aligned_bases = 0;

  friend std::ostream& operator<<(std::ostream& os, const align_count_t& c);
};

struct edge_coverage_t {
  // Read ids that have pair support for this variant branching off from reference
  read_id_set variant_start;
  // Read ids that have pair support for this variant rejoining reference
  read_id_set variant_end;
  // Read ids that have pair support that are in the interior of this assembly.
  read_id_set interior;

  // Read ids that have pair support that counterindicate this variant branching off from reference
  read_id_set reference_start;
  // Read ids that have pair support that counterindicate this variant rejoining reference
  read_id_set reference_end;

  // If sequences are compared between variant and reference, these
  // are the number of shared bases at each end.
  aoffset_t start_common = 0;
  aoffset_t end_common = 0;

  friend std::ostream& operator<<(std::ostream& os, const edge_coverage_t& cov);
};

// Features to pass to machine learning.  These are populated at the
// same time as "report_discovered_assembly".
struct assembly_ml_features {
  int score = 0;
  int refspan = 0;
  int lanch = 0;
  int ranch = 0;
  float refgc = 0;
  float altgc = 0;
  dna_sequence alt_seq;
};

class string_set : public absl::btree_set<std::string> {
 public:
  // Same constructors as btree_set
  using absl::btree_set<std::string>::btree_set;

  // Easy implicit conversion from a vector.
  string_set(const std::vector<std::string>& items) : string_set(items.begin(), items.end()) {}

  // Make phase sets convenient to display
  std::string to_string() const;
  // Like to_string, but doesn't include parents.
  std::string to_string_short() const;
  friend std::ostream& operator<<(std::ostream& os, const string_set& phase_ids) {
    return os << phase_ids.to_string();
  }

  // Make phase sets conveinent to do set operations on.
  string_set& operator+=(const string_set&);
  string_set operator+(const string_set&) const;

  string_set& operator-=(const string_set&);
  string_set operator-(const string_set&) const;

  string_set& operator&=(const string_set&);
  string_set operator&(const string_set&)const;

  bool contains(const std::string& id) const { return count(id); }
};
using phase_set = string_set;

// Stores seqset entries for a path through a seqset as part of a graph.
//
// The base sequence in the path is stored externally.
//
// By itself, a seqset_path will end with seqset->ctx_begin() and
// will start with seqset->find(a prefix of the path sequence).
//
// However, as part of a graph, there will often be more than one end.
struct assemble_options;
// Removes all seqset entries in "rs" that are prefixes of other
// entries in "rs".  This is used when we use push_front_drop on a set
// of ranges, where the prefix is useless and redundant.
void seqset_set_dedup_prefixes(absl::btree_set<seqset_range>& rs);
class seqset_path {
 public:
  seqset_path() = default;
  ~seqset_path();

  void add(aoffset_t offset, seqset_range r);
  void add(aoffset_t offset, const absl::btree_set<seqset_range>& r);

  // Update this path with a new set of ends, and propagate it through.
  void propagate_from_end(const absl::btree_set<seqset_range>& new_ends, dna_slice seq,
                          const assemble_options& opts);

  // Returns starts, e.g. get_marked()[0].  An empty result indicates no data available.
  const absl::btree_set<seqset_range>& starts() const;

  // Return ends, e.g. get_marked()[seq.size()].  An empty result indicates no data available.
  const absl::btree_set<seqset_range>& ends() const;

  const absl::btree_map<aoffset_t, absl::btree_set<seqset_range>>& entries() const {
    return m_entries;
  }

  size_t size() const;

  const absl::btree_set<seqset_range>& mates() const { return m_mates; }

  bool empty() const { return m_entries.empty(); }
  void clear();

  void swap(seqset_path& rhs);

 private:
  static const absl::btree_set<seqset_range> g_empty;

  absl::btree_set<seqset_range> entries_at_offset(aoffset_t offset);

  absl::btree_map<aoffset_t, absl::btree_set<seqset_range>> m_entries;

  absl::btree_set<seqset_range> m_mates;
};

struct assembly {
  assembly();
  assembly(optional_aoffset left_off, optional_aoffset right_off, dna_sequence aseq, size_t asm_id);
  assembly(optional_aoffset left_off, optional_aoffset right_off, dna_sequence aseq);
  assembly(const assembly&);
  assembly(assembly&&);
  assembly& operator=(const assembly&);
  assembly& operator=(assembly&&);
  ~assembly();

  size_t assembly_id = 0;

  // During deduplication, additional assembly ids might be added.
  std::vector<size_t> merged_assembly_ids;

  // Offsets and lengths of various parts of the assembly.  The left
  // anchor spans [left_offset - left_anchor_len, left_offset).  The
  // right anchor spans [right_offset, right_offset +
  // right_anchor_len).
  optional_aoffset left_offset = 0;
  int left_anchor_len = 0;
  optional_aoffset right_offset = 0;
  int right_anchor_len = 0;

  // Contains the bases between [left_offset, right_offset + right_anchor_len)
  dna_sequence seq;

  // Stats from tracing
  unsigned trace_steps = 0;
  unsigned unique_pairs_used = 0;
  unsigned min_overlap = 0;
  unsigned left_anchor_ambiguous_bases = 0;

  // Depth of the most covered assembly at this reference position that isn't this assembly, not
  // covering reference.
  unsigned other_depth = 0;
  unsigned other_pair_depth = 0;
  // Depth of reference opposite this assembly.
  unsigned ref_depth = 0;
  // Number of strands (up to a maximum of
  // assemble_options.max_ploids) that this assembly is present on.
  // For instance, with max_ploids = 2, a value of "0" means this call
  // didn't pass filtering, a value of "1" means "0/1", and a value of
  // "2" means "1/1".
  unsigned strand_count = 0;
  // Scale from 0 to 1
  double genotype_quality = 0;

  // Contains all read ids present supporting this assembly.
  // Generated by the tracing stage.  These reads are all facing left
  // with respect to reference.
  read_id_set rc_read_ids;

  // Interbase coverage (depth); this array is of length seq.size() - 1
  std::vector<int> coverage;

  // Interbase coverage (depth), but only for reads that match pairs.
  std::vector<int> pair_coverage;

  // Read_ids of the reverse complement reads (so that they match the
  // read ids seen using push_front_drop) that provide pairing
  // evidence of this assembly.
  std::vector<uint32_t> left_pair_matches;
  std::vector<uint32_t> right_pair_matches;

  // Resulting assembly score, generated by the pair counting stage.
  acost_t score = 0;

  // True if this assembly matches reference entirely.  Generated by
  // the assembly splitting stages.
  bool matches_reference = false;

  // This contains variants within this assembly that have been
  // aligned.
  std::vector<aligned_var> aligned_variants;

  // True if we should bypass coverage calculations for this assembly.
  bool bypass_coverage = false;

  boost::optional<edge_coverage_t> edge_coverage;
  boost::optional<read_coverage_t> read_coverage;
  boost::optional<read_coverage_t> pair_read_coverage;
  boost::optional<align_count_t> align_count;

  // Minimum value of assemble_options::max_coverage_paths needed
  // to get all read coverages detected for this assembly.
  size_t read_cov_max_paths = 0;

  // Arbitrary data for passing through the pipeline
  boost::any user_data;

  // Arbitrary strings associated with this assembly.  Includes the previous "generated_by" field.
  string_set tags;

  // List of phase IDs associated with this assembly.
  phase_set phase_ids;

  // Subassemblies that are a part of this assembly.

  // TODO(nils): This shared ptr is required so that "assembly" is
  // still copyable.  Is there a different way to do this?  Or can we
  // make assembly_ptr a shared_ptr instead of a unique_ptr?
  std::vector<std::shared_ptr<assembly_ptr>> sub_assemblies;

  // Seqset ranges from push_front_drop(b.complement()) for b in
  // a->seq.rev_comp().  Without surrounding entries,
  // seqset_entries.front() would be seqset->find(a->seq), and
  // seqset_entries.back() would be seqset->ctx_begin().
  //
  // This is a cache.  If seqset_entries isn't empty, it should have
  // exactly seq.size() + 1 entries.  An empty set indicates a lack of
  // data, not a nonexistant path.
  //
  // If seqset_entries[x] and seqset_entries[x + 1] are both non-empty, then
  // seqset_entries[x] must be equal to push_front_drop(seq[x], seqset_entries[x+1].
  seqset_path seqset_entries;
  seqset_path rc_seqset_entries;

  void output_offsets(std::ostream& os) const;
  void output_other_info(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const assembly& as);

  // Sort functions
  using ordering_t = std::function<bool(const assembly&, const assembly&)>;
  static bool left_offset_less_than(const assembly& a, const assembly& b) {
    return a.left_offset < b.left_offset;
  }
  static bool left_anchor_end_less_than(const assembly& a, const assembly& b) {
    return (a.left_offset + a.left_anchor_len) < (b.left_offset + b.left_anchor_len);
  }

  boost::optional<assembly_ml_features> ml_features;
};

enum class sort_order { LEFT_OFFSET_ONLY, OLD_DISCOVER, GRAPH_DISCOVER };
class canon_assembly_order {
 public:
  // TODO(nils): Get rid of the version setting.
  canon_assembly_order(sort_order order = g_default_sort_order) : m_sort_order(order) {}
  // Usable as a sort order.  Guaranteed to put things with earlier left_offsets first.
  bool operator()(const assembly& a, const assembly& b) const;

  // Allow sorting of assemblies by their pointers.
  bool operator()(const assembly_ptr& a, const assembly_ptr& b) const { return (*this)(*a, *b); }
  bool operator()(const assembly* a, const assembly* b) const { return (*this)(*a, *b); }

  static canon_assembly_order from_old_flag(bool old_sort_order) {
    if (old_sort_order && g_default_sort_order == sort_order::OLD_DISCOVER) {
      return canon_assembly_order(sort_order::LEFT_OFFSET_ONLY);
    } else {
      return canon_assembly_order(g_default_sort_order);
    }
  }

  static void set_default_sort_order(sort_order order) { g_default_sort_order = order; }

 private:
  sort_order m_sort_order = sort_order::GRAPH_DISCOVER;
  static sort_order g_default_sort_order;
};

size_t allocate_assembly_id();

struct half_aligned_assembly {
  std::string scaffold_name;
  aoffset_t offset;

  // true if "offset" is the right anchor, false if it's the left anchor.
  bool right_anchor;

  // If right_anchor is false, this half aligned variant is anchored on the left starting at
  // seq.begin().
  // If right_anchor is true, this half aligned variant is anchored on the right starting at
  // seq.end().
  dna_sequence seq;
  size_t assembly_id;
  // Supporting read ids, facing left with respect to reference.
  read_id_set rc_read_ids;
  friend std::ostream& operator<<(std::ostream& os, const half_aligned_assembly& ha);
};

void reverse_assembly_in_place(assembly* a, const readmap* rm, aoffset_t ref_end_pos);

half_aligned_assembly reverse_half_aligned(half_aligned_assembly ha, const readmap* rm,
                                           aoffset_t ref_end_pos);

class scaffold;
struct assemble_stats;

struct assemble_options {
  // Assembly sources
  const ::seqset* seqset = nullptr;
  const ::readmap* readmap = nullptr;
  const reference* ref = nullptr;
  const ref_map* rmap = nullptr;
  const ::variants::scaffold* scaffold = nullptr;
  std::string scaffold_name;

  // Assembly IDs to display additional debug info on.
  std::set<size_t> trace_assembly_ids;

  // Maximum number of rejoins per trace.
  acost_t ambiguous_branch_cost = 1;
  acost_t rejoin_local_cost = 1;
  unsigned max_pairs_per_read = 20;
  acost_t base_cost = 1;
  acost_t increase_max_between_pair_cost = 1;
  acost_t max_cost = 1000 * 1000;
  unsigned max_next_paths = 1024;
  acost_t pairs_used_cost = -100;

  acost_t avg_coverage_score = 300;
  acost_t pair_match_score = 130;
  acost_t min_overlap_score = 100;
  acost_t min_coverage_score = 100;

  acost_t anchor_drop_score = -10000;

  size_t min_anchor_drop_overlap = 15;
  unsigned max_ambiguous_bases = 300;

  acost_t traverse_ref_cost = 100 * 1000;

  acost_t seen_entry_before_cost = 0;  // 500 * 1000;

  bool calculate_coverage = false;

  bool only_trace_forward = false;
  bool trace_dead_ends = true;

  // If nonzero, maximum number of reads that each seqset entry can supply for read coverage.
  unsigned read_cov_max_reads_per_entry = 0;

  // Minimum overlap between reads for tracing.
  size_t min_overlap = 100;

  // Minimum overlap between reads for tracing using the pop_tracer
  size_t min_pop_overlap = 15;

  // Amount of variance allowed in estimated positions when
  // constructing assemblies using the pop tracer.
  size_t pop_tracer_offset_slop = 1000;

  // Number of bases to read ahead to look for local rejoins.  Only
  // deletions smaller than this many bases be detected (but others
  // may possibly be detected as break ends).
  size_t read_ahead_distance = 100000;

  // Size of chunks to split a scaffold into to parallelize.
  size_t scaffold_split_size = 1000 * 1000;

  // Costs, used to trace out rejoin paths:

  // Cost per read that matches an ambiguous location in reference (to
  // cut down on digging through tangles).
  int64_t cost_per_ambiguous_reference = 3000;

  // Number of reads of reference match before calling a break end.
  // cost_per_ref_match * break_end_reads should be more than
  // dead_end_cost so that pairing data is required to call a break
  // end.
  size_t break_end_reads = 200;

  // Cost for each decrease in overlap.
  int64_t decrease_overlap_cost = 1000;

  // Cost for each base of size difference between reference and variant.
  int64_t size_change_cost = 5;

  // Cost for each matching pair found (normally negative to indicate a bonus):
  int64_t pair_match_cost = -20 * 1000;

  // Cost for each non-matching pair found:
  int64_t non_matching_pair_cost = 10 * 000;

  // Cost to rejoin as a dead end:
  int64_t dead_end_cost = 500 * 1000;

  // Maximum number of branches to process between pair verification.
  unsigned max_branches_between_pairs = 2;

  // Maximum number of paths to trace from each reference read.
  size_t max_rejoins = 5;

  // Maximum number of alleles to report.
  size_t max_ploids = 2;

  // Maximum distance between pair verification.
  size_t max_bases_between_pairs = 300;

  // Maximum total search steps per reference read before giving up.
  unsigned max_search_steps = 1000;
  unsigned max_ambiguous_search_steps = 100;

  unsigned initial_search_steps = 100;
  unsigned max_search_steps_per_read = 3;

  // True if forward pairs (is_forward = true) only face inward.
  // False if they only face outwards.
  bool forward_pairs_face_inward = true;

  // Number of bases to look back for a matching pair.
  aoffset_t min_pair_distance = 100;
  aoffset_t max_pair_distance = 1000;

  // If true, ignore pairs associated with any read that matches
  // reference in multiple locations.
  bool ignore_ambiguous_ref_pairs = true;

  // If either of these are true, skip executing the push tracer in
  // that direction.
  bool skip_push_trace_fwd = false;
  bool skip_push_trace_rev = false;

  // If either of these are true, skip executing the pop tracer in
  // that direction.
  bool skip_pop_trace_fwd = false;
  bool skip_pop_trace_rev = false;

  // If true, try to do some simple valid variants filtering and
  // genotyping during the assembly process.
  bool simple_genotype_filter = true;

  // If true, use pop_tracer instead of tracer.
  bool use_pop_tracer = false;

  // If true, use the bidirectional push-pop tracer.
  bool use_bidir_tracer = false;

  // If true, bidir tracer will emit all rejoins found even if paths
  // have been found with a better min overlap.
  bool bidir_tracer_emit_all_rejoins = false;

  // If true, pop trace all reads when doing bidir trace instead of
  // just the ones where overlap decreases.
  bool bidir_pop_all_reads = false;

  // If true, save right partials for all reads instead of just the
  // ones where overlap decreases.
  bool bidir_right_partial_all_reads = false;

  // Number of steps to search on each branch off of reference before giving up.
  size_t bidir_max_branch_steps = 400;

  // Number of extra search steps to allow for each additional pair
  // match found on a branch.
  size_t bidir_branch_steps_per_pair = 100;

  // Don't try to discover more variants after discovering this many
  // alleles in the same place.
  size_t bidir_max_ploids = 3;

  // Absolute minimum anchor length for bidir tracer.  Must not be
  // larger than any minimum overlap.
  size_t bidir_min_anchor_len = 10;

  // Minimum overlap for rejoining to reference without a size change.
  size_t bidir_min_local_ref_overlap = 10;

  // If the portion of the seqset that the popped range encompasses is larger than
  // 1/bidir_max_pop_seqset_portion, don't waste time trying to pop more.
  size_t bidir_max_pop_seqset_portion = 100000;

  // If true, track how much time each branch spends searching and
  // report the slowest branches to search.
  bool bidir_report_slow_branches = false;

  // If true, while tracing, continually validate the invariants in
  // the tracer internal state.  Higher numbers cause more validation,
  // but have more of a performance hit.
  int bidir_validate_trace_state = 0;

  // If true, bidir's pop searching makes a new push search when it
  // finds pair support.  Otherwise, it makes a new pop search.
  bool bidir_pop_makes_push = true;

  // If true, take anchor drop assemblies from the first tracer pass
  // and use those to run a second pass using pop_tracer.
  bool pop_trace_anchor_drop = true;

  // If false, don't start tracing from an ambiguous reference location.
  bool trace_ambiguous_ref = false;

  // Generate a graph of paths when assembling.  This function is
  // called with a "dot" format graph.
  std::function<void(const std::string&)> debug_paths;

  // Factor of the sequence size to use when finding small parts of
  // reference that match.  For instance, with a factor of 4, and the
  // minimum of (sequence length, reference length) is 100, reference
  // matches must be at least 25 bases.
  int ref_align_factor = 4;

  // If we can, always search for reference matches down to this
  // number of bases, even if the calculation using ref_align_factor
  // gives a higher number.
  int max_ref_align_bases = 30;

  // Treat any sequence matching reference as a read for the purposes
  // of calculating overlap.
  bool bidir_treat_ref_as_reads = true;

  // When we're unable to get an anchor on one side of the assembly,
  // this is how many bases of reference to look through compared to
  // the size of the sequence in order to find an anchor that isn't a
  // full read.
  double anchor_drop_size_multiplier = 1.5;

  // If false, reference assemblies are discarded except for when
  // generating genotyping information.
  bool trace_reference_assemblies = false;

  // If true, exclude non-structural variants without pair coverage.
  bool rvg_exclude = false;

  // When genotyping, ignore anything that has less than this
  // portion of depth compared to the best percentage of depth for
  // the corresponding area of reference.
  double min_depth_portion = .23;

  // Throw out calls with less than this much evidence.
  int min_read_depth = 1;

  // Throw out calls with less than this much pairing evidence.
  int min_pair_depth = 0;
  double min_avg_pair_depth = 1;

  // If nonzero, maximum number of simultaneous paths to trace during
  // coverage calculations.
  unsigned max_coverage_paths = 4;

  // If true, penalize depth for sequences that are covered with reads
  // in primarily one direction.  This lets us discount certain types of e.g. Illumina read errors.
  bool penalize_directional_coverage = true;

  // VCF entries with either ALT or REF sequences at least this many
  // bases long will receive additional "SVLEN", "SVTYPE", and "END" fields in
  // the vcf output.
  unsigned vcf_sv_size_threshold = 20;

  // Output AID field for each variant
  bool output_assembly_ids = false;

  // Include ML features in output
  bool output_ml_features = false;

  // When using graph discover, only trace assemblies that have at least one of these tags.
  string_set discover_tags;

  std::function<void(const std::string& /* scaffold name */, double /* seconds */,
                     aoffset_t /* offset */, assemble_stats)>
      report_long_traces_func;

  std::function<void(const assembly& a, bool right_anchor)> report_anchor_drop_func;

  std::function<void(const std::string& /* scaffold name */, aoffset_t /* start */,
                     aoffset_t /* limit */, bool /* rev_comp */, double /* seconds */,
                     assemble_stats)>
      report_chunk_stats_func;

  std::function<void(const assemble_options& options, assembly& a)>
      report_discovered_assemblies_func;

  std::function<void(const assemble_options& options, const assembly& a)>
      report_aligned_assemblies_func;

  std::function<void(const half_aligned_assembly&)> report_half_aligned_func;

  std::function<void(const assemble_options& options, const assembly& a,
                     const std::vector<const assembly*>& better_assemblies)>
      report_genotype_discard_func;

  std::function<void(::variants::discovery::state* st)> report_bidir_initialized_func;

  static const assemble_options& defaults() { return g_defaults; }

 private:
  static assemble_options g_defaults;
};

// Facility for tracing what happens as an assembly goes through the pipeline.
void add_assembly_trace(size_t assembly_id);
// Facility for tracing assembly discovery at a certain offset
void add_offset_trace(aoffset_t offset);
// Reset all traces
void reset_assembly_trace();
extern std::set<size_t> g_trace_assembly_ids;
extern std::set<aoffset_t> g_trace_offsets;
extern bool g_trace_all_assemblies;
extern bool assembly_needs_trace(const assembly& a);
extern bool offset_needs_trace(aoffset_t offset);

// offset_split_pos is relative to left->offset in this verison of
// split_assembly:
std::pair<assembly_ptr /* left */, assembly_ptr /* right */> split_assembly(
    assembly_ptr a, aoffset_t assembly_split_pos, aoffset_t offset_split_pos);
// abs_offset_split_pos is an absolute offset specification of where
// to cut in this version of split_assembly:
std::pair<assembly_ptr /* left */, assembly_ptr /* right */> split_assembly_absoffset(
    assembly_ptr a, aoffset_t assembly_split_pos, optional_aoffset abs_offset_split_pos);

// Pads an assembly with reference bases.
void pad_assembly(assembly* a, aoffset_t new_left_offset, aoffset_t new_right_offset,
                  const assemble_options& options);

// Returns nullptr if there were any conflicts in merging.  Assemblies
// must not be disjoint.
assembly_ptr merge_assemblies(const assembly& a, const assembly& b);

// Checks an assembly for consistency and fails if it is inconsistnet.
void check_assembly(const assembly& a, const std::string from_where);

// Checks an assembly not generated internally, and throws an
// exception if it is inconsistent.
void check_assembly_from_user(const assembly& a);

// An assembly including variants.
std::string dump_assembly_and_vars(const assembly& a);
std::string dump_coverage(const std::vector<int>& cov);

struct assemble_stats : public autostats_base {
  DECLARE_AUTOSTATS(                               //
      assemble_stats, ((COUNTER, ref_reads))       //
      ((COUNTER, ambiguous_ref_reads))             //
      ((COUNTER, step_count))                      //
      ((COUNTER, too_many_steps))                  //
      ((COUNTER, too_many_ambiguous_steps))        //
      ((COUNTER, output_count))                    //
      ((COUNTER, dead_ends))                       //
      ((COUNTER, empty_assemblies))                //
      ((COUNTER, found_pairs))                     //
      ((COUNTER, matched_pairs))                   //
      ((COUNTER, ambiguous_pair_entries))          //
      ((COUNTER, ambiguous_pairs))                 //
      ((COUNTER, unused_next_paths))               //
      ((COUNTER, unused_rejoins))                  //
      ((COUNTER, far_rejoins))                     //
      ((COUNTER, local_rejoins))                   //
      ((COUNTER, loops))                           //
      ((COUNTER, max_branch_cost))                 //
      ((COUNTER, max_branch_cost_between_pairs))   //
      ((COUNTER, too_many_ambiguous))              //
      ((COUNTER, suboptimal_path_prune))           //
      ((COUNTER, ref_assemblies))                  //
      ((COUNTER, too_many_pairs))                  //
      ((COUNTER, too_far_without_pair))            //
      ((COUNTER, extend_ambiguous_rejoin))         //
      ((COUNTER, search_not_fast_enough))          //
      ((COUNTER, next_paths_too_big))              //
      ((COUNTER, too_many_ambiguous_bases))        //
      ((COUNTER, prune_ambiguous_ref))             //
      ((COUNTER, exceeded_branch_limit))           //
                                                   //
      ((COUNTER, rejoin_local_cost))               //
      ((COUNTER, rejoin_far_cost))                 //
      ((COUNTER, dead_end_cost))                   //
      ((COUNTER, ambiguous_branch_cost))           //
      ((COUNTER, decrease_overlap_cost))           //
      ((COUNTER, increase_max_between_pair_cost))  //
      ((COUNTER, base_cost))                       //
      ((COUNTER, pairs_used_cost))                 //
      ((COUNTER, traverse_ref_cost))               //
      ((COUNTER, seen_entry_before_cost))          //
      //
      ((MAX, max_ambiguous_step_count))  //
      ((MAX, max_assembly_len))          //
  );

  std::string as_string() const;
};

class assemble_pipeline_interface {
 public:
  void add(assembly_ptr a);
  // TODO(nils): flush() Should be called before destruction everywhere instead of just some places.
  virtual void flush() {}

  virtual void on_assembly(assembly_ptr a) = 0;
  virtual std::string description() const;

  virtual ~assemble_pipeline_interface() = default;

  // For debugging purposes, verify order on all assemblies received.
  static void global_set_verify_order(bool verify_order) { g_verify_order = verify_order; }

 protected:
  void set_expected_order(assembly::ordering_t ordering) { m_expected_order = ordering; }
  static bool g_verify_order;

 private:
  assembly::ordering_t m_expected_order;
  assembly_ptr m_last_assembly;
};

using pipeline_step_t = std::unique_ptr<assemble_pipeline_interface>;

class sorted_output_pipeline_step : public assemble_pipeline_interface {
 protected:
  sorted_output_pipeline_step(pipeline_step_t output, bool old_sort_order = false)
      : m_output(std::move(output)),
        m_output_queue(canon_assembly_order::from_old_flag(old_sort_order)) {}
  ~sorted_output_pipeline_step();

  void flush_sorted_to(aoffset_t offset);
  void flush_sorted();
  void sort_and_output(assembly_ptr a);

  // Track currently in-progress left offsets so our output is sorted.
  void track_left_offset(aoffset_t offset);
  void untrack_left_offset(aoffset_t offset);
  aoffset_t sort_flush_point() const { return m_flush_point; }

 protected:
  std::string sorted_output_stats(boost::optional<aoffset_t> relative_to = boost::none) const;

 private:
  aoffset_t m_flush_point = std::numeric_limits<aoffset_t>::min();
  pipeline_step_t m_output;

  // Have to have a non-move-only pointer to assembly as the key since otherwise there's
  // no way to extract elements from the queue.
  std::multimap<assembly*, assembly_ptr, canon_assembly_order> m_output_queue;
  std::multiset<aoffset_t> m_left_offsets;
};

class assemble_lambda_output : public assemble_pipeline_interface {
 public:
  assemble_lambda_output(const std::function<void(assembly_ptr)>& output_f,
                         const std::string& description);
  assemble_lambda_output() = delete;

  void on_assembly(assembly_ptr a) override;
  std::string description() const override { return m_description; }

 private:
  std::function<void(assembly_ptr)> m_output_f;
  std::string m_description;
};

class assemble_lambda_copy : public assemble_pipeline_interface {
 public:
  assemble_lambda_copy(const std::function<void(const assembly&)>& copy_f, pipeline_step_t output,
                       const std::string& description);
  assemble_lambda_copy() = delete;

  void on_assembly(assembly_ptr a) override;
  std::string description() const override { return m_description; }

 private:
  std::function<void(const assembly&)> m_copy_f;
  pipeline_step_t m_output;
  std::string m_description;
};

class pipeline_interface {
 public:
  virtual ~pipeline_interface() = default;
  virtual pipeline_step_t make_parallel_input() = 0;

 protected:
  pipeline_interface() = default;
};

class scaffold_pipeline_interface {
 public:
  virtual ~scaffold_pipeline_interface() = default;

  virtual std::unique_ptr<pipeline_interface> pipeline_for_scaffold(
      const assemble_options& options, const std::string& scaffold_name) = 0;

 protected:
  scaffold_pipeline_interface() = default;
};

// Hash function that doesn't vary between runs; std::hash makes no guarantee for this.
struct unsalted_hash {
  size_t operator()(size_t in) const { return in; }
};

}  // namespace variants

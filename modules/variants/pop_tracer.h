#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

class pop_tracer {
 public:
  static const char k_pop_tracer_name[];

  pop_tracer(const assemble_options& options);
  pop_tracer() = delete;
  pop_tracer(const pop_tracer&) = delete;
  ~pop_tracer();

  // Adds the given section of reference as anchors.
  void add_reference(aoffset_t start, aoffset_t limit);

  // Adds the given read as a potential unanchored read to be used for
  // assemblies within the given offset range.  "description" is only
  // used for debug tracing.
  void add_read(uint32_t read_id, aoffset_t start_offset, aoffset_t limit_offset);

  void add_anchor_drop(const assembly& ha, bool right_anchor);

  void assemble(assemble_pipeline_interface* output);

  static void add_debug_read(uint32_t read_id);
  static void add_debug_seqset_entry(const seqset_range& r);
  static void clear_debug_reads();

 private:
  // Step 1: Add all reference seqset entries, with their original ref
  // location.  These entries will have both a "left" and "right"
  // offset.
  //
  // Step 2: Get seqset entries for all mates, along with their range
  // of possible locations.
  //
  // Step 3: For all of the entries in 1 and 2, save in a lookup table "fronts".
  //
  // Step 4: For all the entries in 1 and 2, save in a table "poppers".
  //
  // Step 5: For each entry in "poppers":
  //   5a: If the popped entry is less than min_overlap bases, discard this new popper.
  //   5b: If the popped entry is a prefix of any entry in the lookup table "fronts", then:
  //    * Consume the matched entry from "fronts", removing it from
  //     "fronts".
  //    * The new entry has a left offset from "fronts" if present,
  //      and a right offset from the "popper".
  //    * If this new entry has both a left offset and a right offset, an assembly can be emitted,
  //      subsequently discarding this popper.
  //    * Save this popped entry in "new poppers".
  //
  // Step 6: If any entries made it to "new poppers", move it back to "poppers" and repeat step 5.

  struct entry {
    // If present, the known left offset of this entry's anchor in
    // reference.  This is at the head of orig_r.
    boost::optional<aoffset_t> left_offset;
    boost::optional<aoffset_t> right_offset;

    // Reference range in which this entry can be placed.
    aoffset_t start_limit, end_limit;

    // Original range, facing towards the right
    seqset_range orig_r;

    // Popped range, facing towards the right
    seqset_range popped_r;

    // Sequence of bases which has been popped.  Invariant:
    // popped_r.push_front_drop(seq) == orig_r
    dna_sequence seq;

    std::vector<uint32_t> seen_read_ids;

    bool matches_reference = false;
    bool trace_this = false;
  };
  friend std::ostream& operator<<(std::ostream& os, const entry& p);

  struct seqset_range_comparer {
    bool operator()(const seqset_range& lhs, const seqset_range& rhs) const {
      if (lhs.begin() != rhs.begin()) {
        return lhs.begin() < rhs.begin();
      }
      return lhs.size() < rhs.size();
    }
  };
  struct match_sorter;
  struct popper_sorter;

  void add_popper(std::shared_ptr<entry> p);
  void add_reference_at(aoffset_t offset, const seqset_range& r);
  void pop_pass();
  void match_and_output_pass(assemble_pipeline_interface* output);
  void output_assembly(assemble_pipeline_interface* output, const entry& p);
  void add_entry_reads(entry& p, const seqset_range& r);

  static bool range_needs_trace(const seqset_range& r);
  static bool entry_needs_trace(const entry& p);

  using fronts_t = std::multimap<seqset_range, std::shared_ptr<entry>, seqset_range_comparer>;
  fronts_t m_fronts;
  std::vector<std::shared_ptr<entry>> m_poppers;

  assemble_options m_options;
};

}  // namespace variants

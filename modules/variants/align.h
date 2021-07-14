#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

class aligner : public sorted_output_pipeline_step {
 public:
  aligner(const assemble_options& options, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output), true /* old sort order */),
        m_options(options),
        m_scaffold(options.scaffold) {
    CHECK(m_scaffold);
    set_expected_order(assembly::left_offset_less_than);
  }

  void on_assembly(assembly_ptr a) override;

  // Method to use to score matches.
  //
  // ANCHORED_TO_LEFT means that the right anchor has been dropped, so
  // try to get the match as close to the same position as possible
  // when the left edges have been liened up.
  //
  // ANCHORED_TO_RIGHT is the same, but for left anchor drops.
  //
  // ANCHORED_TO_BOTH means to get the match as close to the center as
  // possible.
  enum class anchor_type_t {
    ANCHORED_TO_LEFT,
    ANCHORED_TO_RIGHT,
    ANCHORED_TO_BOTH
  };

  // Return true if a match is found.
  static bool find_biggest_match(const assemble_options& options, dna_slice seq,
                                 const scaffold& s, int* match_len,
                                 aoffset_t* seq_match_start,
                                 aoffset_t* scaffold_match_start,
                                 aoffset_t* min_match_size,
                                 anchor_type_t anchor);
  static bool find_biggest_match_with_ends(const assemble_options& options,
                                           dna_slice seq, const scaffold& s,
                                           int* match_len,
                                           aoffset_t* seq_match_start,
                                           aoffset_t* scaffold_match_start);

  static bool find_end_matches(const assemble_options& options, dna_slice seq,
                               const scaffold& s, int* match_len,
                               aoffset_t* seq_match_start,
                               aoffset_t* scaffold_match_start,
                               aoffset_t max_match_size);

  // Returns true if a match of the given size has found.
  static bool find_match(dna_slice seq, const scaffold& s, int match_len,
                         aoffset_t* seq_match_start,
                         aoffset_t* scaffold_match_start, anchor_type_t anchor);

 private:
  struct work_item {
    aligned_var v;
    scaffold s;
  };

  void process(assembly* a, work_item w);

  assemble_options m_options;
  const scaffold* m_scaffold = nullptr;

  std::vector<work_item> m_work;
};

class align_splitter : public sorted_output_pipeline_step {
 public:
  align_splitter(pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output), true /* old sort order */) {
    set_expected_order(assembly::left_offset_less_than);
  }
  void on_assembly(assembly_ptr a) override;

 private:
  void set_matches_reference(assembly& a);
};

std::ostream& operator<<(std::ostream& os,
                         const aligner::anchor_type_t& anchor);


}  // namespace variants

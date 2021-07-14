#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

// * Trims all sections of incoming non-reference assemblies that
//   match reference, and shrinks the assemblies.
//
// * Slices up reference assemblies so that there aren't any non-reference anchors that
//   output in the middle of reference assemblies.
class graph_trim_ref : public sorted_output_pipeline_step {
 public:
  graph_trim_ref(const assemble_options& options, pipeline_step_t output);
  graph_trim_ref() = delete;
  ~graph_trim_ref() override;

  void on_assembly(assembly_ptr a) override;
  void flush() override;

 private:
  void advance_to(aoffset_t pos);
  void split_and_output_ref(assembly_ptr a);

  assemble_options m_options;
  const scaffold* m_scaffold = nullptr;

  // Maximum number of bases to backtrack when trimming reference
  // bases on the right side of a half-anchored assembly.
  // TODO(nils): Should this be configurable?
  aoffset_t m_max_backtrack = 300;

  // Reference assemblies
  std::deque<assembly_ptr> m_ref_asms;

  // Offsets where we need to make sure reference assemblies stop.
  absl::btree_set<aoffset_t> m_ref_stops;
};

}  // namespace variants

#pragma once

#include "modules/bio_base/seqset.h"
#include "modules/variants/assemble.h"

namespace variants {

// add_ref does not change assemblies that pass through it, but adds
// reference assemblies covering all sections of reference that the
// assemblies cover.  These can be used to compare coverage between
// variants and reference.
//
// Note that there may be reference offsets that incoming assemblies
// cover that outgoing reference assemblies do not in the case of
// missing ("NNNNNN...") sections of reference.  In this case,
// reference assemblies covering those missing regions will not be
// emitted.

class add_ref : public sorted_output_pipeline_step {
 public:
  // Add pad_size bases of reference around every assembly.  When
  // generating coverage, this can be used to make sure we trace
  // enough to get pairs and full reads on each side.
  //
  // If whole_ref is true, add whole reference even if no input assemblies are present.
  add_ref(const assemble_options& options, aoffset_t pad_size, bool whole_ref, int max_len,
          pipeline_step_t output);
  add_ref() = delete;
  ~add_ref() { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  void flush();
  void advance_to(aoffset_t offset);
  void output_ref(aoffset_t start, aoffset_t limit);
  void output_ref_part(aoffset_t start, aoffset_t limit);

  // Offsets where we need to start/stop reference assemblies.
  std::set<aoffset_t> m_edge_offsets;

  // Current offset where we're outputting reference assemblies to.
  aoffset_t m_cur_offset = 0;

  aoffset_t m_padded_right_offset = 0;

  aoffset_t m_pad_size = 0;
  aoffset_t m_max_len = 0;

  assemble_options m_options;

  static const char k_add_ref_name[];
};

}  // namespace variants

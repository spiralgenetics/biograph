#pragma once

#include "modules/variants/assemble.h"

namespace variants {

// Genotypes input assemblies based on the min depth of the variants.
class genotyper : public sorted_output_pipeline_step {
 public:
  genotyper(const assemble_options& options, pipeline_step_t output);
  virtual ~genotyper() { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  struct entry {
    assembly_ptr a;

    // True if this entry is an active allele.
    bool active = false;
    std::vector<aligned_var>::iterator vit;

    aoffset_t seq_offset = 0;
    aoffset_t ref_offset = 0;

    // This is where we split if we're deactivating at ref_offset - 1.
    aoffset_t deactivate_seq_offset = 0;
    aoffset_t deactivate_ref_offset = 0;

    // Depth of variant pointed to by "vit":
    int variant_depth = 0;

    // Current depth to the left of the base at m_process_offset.
    int cur_depth = 0;

    // True if we're in the middle of a variant and can't activate.
    bool in_variant = false;

  };
  friend std::ostream& operator<<(std::ostream& os, const entry& e);
  using entry_ptr = std::unique_ptr<entry>;

  void init_entry(entry& e);
  void advance();
  void advance_base();
  void advance_entry(entry& entry);
  void calc_depth(entry& entry);
  void flush();

  void add_assembly(assembly_ptr a);

  bool check_entry_done(entry& e);

  void activate(entry& e);
  void deactivate(entry& e);

  bool is_degenerate(const assembly& a);
  void report_discard(assembly_ptr a);

  std::vector<std::unique_ptr<entry>> m_entries;

  aoffset_t m_process_offset = std::numeric_limits<aoffset_t>::min();
  aoffset_t m_intake_offset = std::numeric_limits<aoffset_t>::min();

  assemble_options m_options;
};

}  // namespace variants

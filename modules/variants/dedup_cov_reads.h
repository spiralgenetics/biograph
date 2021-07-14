#pragma once

#include "modules/variants/assemble.h"

namespace variants {

// Deduplicates edge coverage reads so that each read is only present
// once.  This allows us to better genotype repetitive regions.
//
// This must come after pair_edge_cov in the pipeline.

class dedup_cov_reads : public sorted_output_pipeline_step {
 public:
  dedup_cov_reads(const assemble_options& opts, pipeline_step_t output);
  ~dedup_cov_reads() override;

  void on_assembly(assembly_ptr a) override;

 private:
  void advance_to(aoffset_t offset);
  void track_reads(assembly* a, bool track /* true = track, false = untrack */);
  void track_ref(assembly* a, read_id_set* reads, bool track );
  void track_var(assembly* a, read_id_set* reads, bool track );
  void remove_read_if_seen_in_var(assembly* because_of, uint32_t read_id);
  void flush();

  aoffset_t m_cur_offset = 0;
  assemble_options m_options;

  struct var_seen_t {
    read_id_set* collection = nullptr;
    assembly* a = nullptr;
  };

  std::unordered_map<uint32_t, unsigned /* seen count */, unsalted_hash> m_seen_ref_reads;

  std::multimap<aoffset_t /* right offset */, assembly_ptr> m_active;

  std::unordered_multimap<uint32_t /* read id */, var_seen_t, unsalted_hash> m_seen_var_reads;
};

}  // namespace variants

#pragma once

// Finds and removes reads that have multiple alignments between adjacent assemblies.

#include "modules/variants/apply_block.h"
#include "modules/variants/apply_edges.h"
#include "modules/variants/assemble.h"

namespace variants {

class filter_dup_align : public apply_block_step {
 public:
  using sort_func_t = std::function<std::vector<assembly_ptr>(std::vector<assembly_ptr>)>;
  filter_dup_align(const sort_func_t& sort_func, pipeline_step_t output);
  ~filter_dup_align() override;
  void on_block(aoffset_t left_offset, aoffset_t right_offset,
                const std::vector<assembly_ptr>& block) override;

 private:
  struct trace_state;
  struct search_entry;
  class block_tracer;

  sort_func_t m_sort_func;
};

}  // namespace variants

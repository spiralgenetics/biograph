#pragma once

#include "modules/variants/assemble.h"

namespace variants {

class apply_block_step : public sorted_output_pipeline_step {
 public:
  apply_block_step(pipeline_step_t output);
  ~apply_block_step();

  void on_assembly(assembly_ptr a) override;

 protected:
  virtual void on_block(aoffset_t left_offset, aoffset_t right_offset,
                        const std::vector<assembly_ptr>& block) = 0;

  // Must be called from the destructor of subclasses before they destroy themselves.
  void flush();

 private:
  void advance_to(aoffset_t offset);
  void process_current();
  void flush_block();

  aoffset_t m_block_start = std::numeric_limits<aoffset_t>::min();
  aoffset_t m_block_end = std::numeric_limits<aoffset_t>::min();
  aoffset_t m_cur_offset = std::numeric_limits<aoffset_t>::min();

  // True if multiple assemblies have been seen ending on m_block_end.
  bool m_block_end_multiple = false;

  // Assemblies starting at m_cur_offset
  std::vector<assembly_ptr> m_current;

  // Assemblies in the current block waiting to be processed.
  std::vector<assembly_ptr> m_block;
};

using apply_block_func_t = std::function<void(
    aoffset_t /* left_offset */, aoffset_t /* right_offset */, const std::vector<assembly_ptr>&)>;

class apply_block_lambda_step : public apply_block_step {
 public:
  apply_block_lambda_step(pipeline_step_t output, const apply_block_func_t& on_block);
  ~apply_block_lambda_step() override;
  void on_block(aoffset_t left_offset, aoffset_t right_offset,
                const std::vector<assembly_ptr>& block) override;

 private:
  const apply_block_func_t& m_on_block;
};

}  // namespace variants

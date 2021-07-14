#pragma once

#include "modules/variants/assemble.h"

namespace variants {

class apply_edges_step : public sorted_output_pipeline_step {
 public:
  apply_edges_step(pipeline_step_t output);
  ~apply_edges_step();

  void on_assembly(assembly_ptr a) override;

  void flush() override;

 protected:
  virtual void on_assembly_edges(optional_aoffset reference_pos,
                                 const std::vector<assembly_ptr>& left_edges,
                                 const std::vector<assembly_ptr>& inserts,
                                 const std::vector<assembly_ptr>& right_edges) = 0;
  virtual void on_advance(aoffset_t new_cur_offset) {}

 private:
  void advance_to(aoffset_t offset);
  void advance_towards(aoffset_t offset);
  void flush_active_to_here();

  aoffset_t m_cur_offset = std::numeric_limits<aoffset_t>::min();

  std::vector<assembly_ptr> m_cur_inserts;
  std::vector<assembly_ptr> m_cur_non_inserts;
  std::vector<assembly_ptr> m_cur_left_unanchored;

  std::multimap<aoffset_t /* right offset */, assembly_ptr> m_active;
};

using apply_edges_func_t =
    std::function<void(optional_aoffset reference_pos, const std::vector<assembly_ptr>&,
                       const std::vector<assembly_ptr>&, const std::vector<assembly_ptr>&)>;

class apply_edges_lambda_step : public apply_edges_step {
 public:
  apply_edges_lambda_step(pipeline_step_t output, const apply_edges_func_t& on_edges);
  ~apply_edges_lambda_step() override;
  void on_assembly_edges(optional_aoffset reference_pos,
                         const std::vector<assembly_ptr>& left_edges,
                         const std::vector<assembly_ptr>& inserts,
                         const std::vector<assembly_ptr>& right_edges) override;

 private:
  const apply_edges_func_t& m_on_edges;
};

extern void apply_edges_to_block(std::vector<assembly_ptr>& block,
                                 const apply_edges_func_t& on_block);

}  // namespace variants

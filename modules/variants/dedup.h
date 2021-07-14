#pragma once

#include "modules/variants/assemble.h"

namespace variants {

// Deduplicates assemblies, combining their counters.  Input should be
// sorted by left offset.  Sorting order is preserved.
class deduper : public assemble_pipeline_interface {
 public:
  deduper(pipeline_step_t output) : m_output(std::move(output)) {
    set_expected_order(assembly::left_offset_less_than);
  }
  deduper() = delete;
  ~deduper() { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  assembly_ptr combine(assembly_ptr a, assembly_ptr b);
  void advance_to(aoffset_t offset);
  void flush();

  pipeline_step_t m_output;

  std::multimap<aoffset_t /* left offset */, assembly_ptr> m_queued;
};

class exact_deduper : public assemble_pipeline_interface {
 public:
  exact_deduper(pipeline_step_t output) : m_output(std::move(output)) {
    set_expected_order(assembly::left_offset_less_than);
  }
  exact_deduper() = delete;
  ~exact_deduper() override { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  void advance_to(aoffset_t offset);
  void flush();

  pipeline_step_t m_output;

  std::multimap<aoffset_t /* left offset */, assembly_ptr> m_queued;
};

}  // namespace variants

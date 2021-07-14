#pragma once

#include "modules/variants/assemble.h"

namespace variants {

class sorter : public assemble_pipeline_interface {
 public:
  sorter(const assembly::ordering_t& compare_f, pipeline_step_t output)
      : m_compare_f(compare_f), m_output(std::move(output)) {}
  sorter() = delete;
  ~sorter() {
    flush();
    CHECK(m_queued.empty());
  }

  void on_assembly(assembly_ptr a) override { m_queued.emplace_back(std::move(a)); }

 private:
  void flush();

  assembly::ordering_t m_compare_f;
  pipeline_step_t m_output;

  // Owned pointers, since parallel sort deals poorly with unique ptrs.
  std::vector<assembly_ptr> m_queued;
};

}  // namespace variants

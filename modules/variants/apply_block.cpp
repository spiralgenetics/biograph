#include "modules/variants/apply_block.h"

namespace variants {

namespace {

constexpr bool k_dbg = false;

}  // namespace

apply_block_step::apply_block_step(pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)) {}

apply_block_step::~apply_block_step() {
  CHECK_EQ(m_cur_offset, std::numeric_limits<aoffset_t>::max())
      << "Subclasses of apply_block_step must call flush() in their destructors.";
  CHECK(m_block.empty())
      << "Subclasses of apply_block_step must call flush() in their destructors.";
}

void apply_block_step::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  process_current();
}

void apply_block_step::on_assembly(assembly_ptr a) {
  advance_to(a->left_offset);
  m_current.push_back(std::move(a));
}

void apply_block_step::advance_to(aoffset_t cur) {
  if (cur > m_cur_offset) {
    process_current();
    m_cur_offset = cur;
  } else {
    CHECK_EQ(cur, m_cur_offset);
  }
}

void apply_block_step::process_current() {
  if (k_dbg) {
    std::cerr << "Processing " << m_current.size() << " assemblies at " << m_cur_offset << "\n";
  }
  if (m_cur_offset > m_block_end ||
      (m_cur_offset == m_block_end && !m_block_end_multiple && m_current.size() <= 1)) {
    flush_block();
  }

  for (auto& a : m_current) {
    if (a->right_offset > m_block_end) {
      if (k_dbg) {
        std::cerr << "Extending block end to " << m_block_end << "\n";
      }
      m_block_end = a->right_offset;
      m_block_end_multiple = false;
    } else if (a->right_offset == m_block_end) {
      m_block_end_multiple = true;
    }
    if (m_block.empty()) {
      if (k_dbg) {
        std::cerr << "Starting new block at " << a->left_offset << "\n";
      }
      m_block_start = a->left_offset;
    }
    m_block.emplace_back(std::move(a));
  }
  m_current.clear();
}

void apply_block_step::flush_block() {
  if (k_dbg) {
    std::cerr << "Flushing block of [" << m_block_start << ", " << m_block_end << ") size "
              << m_block.size() << " at " << m_cur_offset << "\n";
  }
  if (m_block.empty()) {
    if (k_dbg) {
      std::cerr << "Flushing sorted to cur " << m_cur_offset << "\n";
    }
    flush_sorted_to(m_cur_offset);
    return;
  }
  on_block(m_block_start, m_block_end, m_block);
  for (auto& a : m_block) {
    sort_and_output(std::move(a));
  }
  m_block.clear();
  if (k_dbg) {
    std::cerr << "Flushing sorted after block to " << m_block_end << "\n";
  }
  flush_sorted_to(m_block_end);
}

apply_block_lambda_step::apply_block_lambda_step(pipeline_step_t output,
                                                 const apply_block_func_t& on_block)
    : apply_block_step(std::move(output)), m_on_block(on_block) {}

apply_block_lambda_step::~apply_block_lambda_step() { flush(); }

void apply_block_lambda_step::on_block(aoffset_t left_offset, aoffset_t right_offset,
                                       const std::vector<assembly_ptr>& block) {
  return m_on_block(left_offset, right_offset, block);
}

}  // namespace variants

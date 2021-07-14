#include "modules/variants/apply_edges.h"

namespace variants {

apply_edges_step::apply_edges_step(pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)) {}

apply_edges_step::~apply_edges_step() {
  CHECK_EQ(m_cur_offset, std::numeric_limits<aoffset_t>::max())
      << "Subclasses of apply_edges_step must call flush() in their destructors.";
  CHECK(m_cur_inserts.empty());
  CHECK(m_cur_non_inserts.empty());
  CHECK(m_active.empty());
  CHECK(m_cur_left_unanchored.empty());
}

void apply_edges_step::on_assembly(assembly_ptr a) {
  aoffset_t pos = min(a->left_offset, a->right_offset);
  track_left_offset(pos);
  advance_to(pos);

  if (a->left_offset == a->right_offset) {
    m_cur_inserts.emplace_back(std::move(a));
  } else {
    if (a->left_offset) {
      m_cur_non_inserts.emplace_back(std::move(a));
    } else {
      m_cur_left_unanchored.emplace_back(std::move(a));
    }
  }
}

void apply_edges_step::advance_to(aoffset_t offset) {
  while (m_cur_offset < offset) {
    advance_towards(offset);
    flush_sorted_to(m_cur_offset);
  }
  CHECK_EQ(m_cur_offset, offset);
}

void apply_edges_step::advance_towards(aoffset_t target_offset) {
  CHECK_GT(target_offset, m_cur_offset);
  flush_active_to_here();
  if (!m_active.empty() && target_offset > m_active.begin()->first) {
    target_offset = m_active.begin()->first;
  }

  CHECK(m_cur_inserts.empty());
  CHECK(m_cur_non_inserts.empty());
  CHECK(m_cur_left_unanchored.empty());
  m_cur_offset = target_offset;

  on_advance(m_cur_offset);
}

void apply_edges_step::flush_active_to_here() {
  if (!m_active.empty()) {
    CHECK_GE(m_active.begin()->first, m_cur_offset);
  }

  std::vector<assembly_ptr> cur_left_assemblies;
  if (!m_cur_left_unanchored.empty()) {
    on_assembly_edges(optional_aoffset::none, {}, {}, m_cur_left_unanchored);
    cur_left_assemblies = std::move(m_cur_left_unanchored);
    m_cur_left_unanchored.clear();
  }

  while (!m_active.empty() && m_active.begin()->first == m_cur_offset) {
    cur_left_assemblies.emplace_back(std::move(m_active.begin()->second));
    m_active.erase(m_active.begin());
  }

  if (!cur_left_assemblies.empty() || !m_cur_inserts.empty() || !m_cur_non_inserts.empty()) {
    on_assembly_edges(m_cur_offset, cur_left_assemblies, m_cur_inserts, m_cur_non_inserts);
  }

  for (auto& a : cur_left_assemblies) {
    CHECK_EQ(a->right_offset, m_cur_offset)
        << "Right offset changed during edge processing: " << *a;
    aoffset_t left_offset = min(a->left_offset, a->right_offset);
    untrack_left_offset(left_offset);
    sort_and_output(std::move(a));
  }
  cur_left_assemblies.clear();

  for (auto& a : m_cur_inserts) {
    CHECK_EQ(a->left_offset, m_cur_offset)
        << "Left offset of insert changed during edge processing: " << *a;
    CHECK_EQ(a->right_offset, m_cur_offset)
        << "Right offset of insert changed during edge processing: " << *a;
    aoffset_t left_offset = a->left_offset;
    untrack_left_offset(left_offset);
    sort_and_output(std::move(a));
  }
  m_cur_inserts.clear();

  std::vector<assembly_ptr> cur_right_unanchored;
  for (auto& a : m_cur_non_inserts) {
    CHECK_EQ(a->left_offset, m_cur_offset) << "Left offset changed during edge processing: " << *a;
    if (a->right_offset) {
      auto right_offset = a->right_offset;
      m_active.emplace(right_offset, std::move(a));
    } else {
      cur_right_unanchored.emplace_back(std::move(a));
    }
  }
  m_cur_non_inserts.clear();

  if (!cur_right_unanchored.empty()) {
    on_assembly_edges(optional_aoffset::none, cur_right_unanchored, {}, {});
    for (auto& a : cur_right_unanchored) {
      aoffset_t left_offset = a->left_offset;
      untrack_left_offset(left_offset);
      sort_and_output(std::move(a));
    }
    cur_right_unanchored.clear();
  }

  if (!m_active.empty()) {
    CHECK_GE(m_active.begin()->first, m_cur_offset);
  }
}

void apply_edges_step::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
  flush_sorted();
}

apply_edges_lambda_step::apply_edges_lambda_step(pipeline_step_t output,
                                                 const apply_edges_func_t& on_edges)
    : apply_edges_step(std::move(output)), m_on_edges(on_edges) {}

apply_edges_lambda_step::~apply_edges_lambda_step() { flush(); }
void apply_edges_lambda_step::on_assembly_edges(optional_aoffset reference_pos,
                                                const std::vector<assembly_ptr>& left_edges,
                                                const std::vector<assembly_ptr>& inserts,
                                                const std::vector<assembly_ptr>& right_edges) {
  m_on_edges(reference_pos, left_edges, inserts, right_edges);
}

void apply_edges_to_block(std::vector<assembly_ptr>& input, const apply_edges_func_t& on_block) {
  std::vector<assembly_ptr> output;
  output.reserve(input.size());

  {
    auto save_output = make_unique<assemble_lambda_output>(
        [&output](assembly_ptr a) { output.emplace_back(std::move(a)); },
        "apply_edges_to_block output");
    apply_edges_lambda_step apply_lambda(std::move(save_output), on_block);
    for (auto& a : input) {
      apply_lambda.add(std::move(a));
    }
  }
  input = std::move(output);
}

}  // namespace variants

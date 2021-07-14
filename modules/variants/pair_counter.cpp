#include "modules/variants/pair_counter.h"
#include "modules/bio_base/readmap.h"

namespace variants {

constexpr int k_pair_counter_debug = 0;

pair_counter::pair_counter(const assemble_options& options, pipeline_step_t output)
    : m_options(options), m_output(std::move(output)) {}

void pair_counter::on_assembly(assembly_ptr a) {
  if (k_pair_counter_debug) {
    std::cout << "Adding unscored assembly: " << *a << "\n";
  }
  if (!m_right.empty()) {
    CHECK_LE(m_right.rbegin()->first, a->left_offset) << "Left offsets should be non-decreasing.";
  }

  while (advance_some(a->left_offset - m_options.max_pair_distance)) {
  }

  for (const auto& read_id : a->rc_read_ids) {
    m_right_read_ids.insert(read_id);
  }
  auto left_offset = a->left_offset;
  m_right.emplace(std::make_pair(left_offset, std::move(a)));
}

void pair_counter::calc_score(assembly& a) const {
  int min_coverage = 0;
  acost_t avg_coverage = 0;
  if (!a.coverage.empty()) {
    auto first_coverage = a.coverage.begin();
    while (first_coverage + 1 != a.coverage.end() && *first_coverage <= *(first_coverage + 1)) {
      ++first_coverage;
    }
    auto last_coverage = a.coverage.end();
    --last_coverage;
    while (last_coverage != first_coverage && *last_coverage <= *(last_coverage - 1)) {
      --last_coverage;
    }

    ++last_coverage;
    int coverage_count = 0;
    int min_coverage = std::numeric_limits<int>::max();
    for (auto it = first_coverage; it != last_coverage; ++it) {
      if (min_coverage > *it) {
        min_coverage = *it;
      }
      avg_coverage += *it;
      coverage_count++;
    }
    avg_coverage /= coverage_count;
  }
  a.score += min_coverage * m_options.min_coverage_score;
  a.score += avg_coverage * m_options.avg_coverage_score;
  a.score += a.min_overlap * m_options.min_overlap_score;
  a.score +=
      (a.left_pair_matches.size() + a.right_pair_matches.size()) * m_options.pair_match_score;
  if (k_pair_counter_debug) {
    std::cout << "Calculated score " << a.score << " for " << a << " including "
              << a.left_pair_matches.size() << " left and " << a.right_pair_matches.size()
              << " right pair matches\n";
  }
}

void pair_counter::flush() {
  while (advance_some(std::numeric_limits<aoffset_t>::max())) {
  }
}

bool pair_counter::advance_some(aoffset_t target_offset) {
  if (m_cur_offset >= target_offset) {
    return false;
  }

  aoffset_t advance_amount = target_offset - m_cur_offset;
  if (k_pair_counter_debug) {
    std::cout << "Advancing " << advance_amount << " from " << m_cur_offset << " to "
              << target_offset << "; " << m_left.size() << " in left, " << m_active.size()
              << " active, " << m_right.size() << " in right\n";

    if (k_pair_counter_debug > 1) {
      std::cout << m_left.size() << " lefts:\n";
      for (const auto& l : m_left) {
        std::cout << " " << l.first << ": " << l.second.size() << " reads\n";
      }
      std::cout << m_active.size() << " active:\n";
      for (const auto& a : m_active) {
        std::cout << " " << a.first << ": " << *a.second << "\n";
      }
      std::cout << m_right.size() << " right:\n";
      for (const auto& r : m_right) {
        std::cout << " " << r.first << ": " << *r.second << "\n";
      }
      std::cout << "\n";
    }
  }

  if (!m_left.empty()) {
    aoffset_t left_advance_limit =
        (m_left.begin()->first + m_options.max_pair_distance) - m_cur_offset + 1;
    advance_amount = std::min(advance_amount, left_advance_limit);
    if (k_pair_counter_debug) {
      std::cout << "Left advance limit: " << left_advance_limit
                << " (first left = " << m_left.begin()->first << ")\n";
    }
  }

  if (!m_active.empty()) {
    aoffset_t active_advance_limit = m_active.begin()->first - m_cur_offset + 1;
    advance_amount = std::min(advance_amount, active_advance_limit);
    if (k_pair_counter_debug) {
      std::cout << "Active advance limit: " << active_advance_limit
                << " (first active = " << m_active.begin()->first << ")\n";
    }
  }

  if (!m_right.empty()) {
    aoffset_t right_advance_limit = m_right.begin()->first - m_cur_offset + 1;
    advance_amount = std::min(advance_amount, right_advance_limit);
    if (k_pair_counter_debug) {
      std::cout << "Right advance limit: " << right_advance_limit
                << " (first right = " << m_right.begin()->first << ")\n";
    }
  }

  CHECK_GT(advance_amount, 0) << "Advance didn't do anything!";
  m_cur_offset += advance_amount;

  advance_right();
  advance_active();
  advance_left();

  return true;
}

void pair_counter::advance_active() {
  while (!m_active.empty() && m_active.begin()->first < m_cur_offset) {
    auto next_active = m_active.begin();
    assembly_ptr a = std::move(next_active->second);
    m_active.erase(next_active);

    a->right_pair_matches =
        find_pair_matches(*a, m_right_read_ids, !m_options.forward_pairs_face_inward);

    read_id_set read_ids = std::move(a->rc_read_ids);
    a->rc_read_ids.clear();

    for (const auto& read_id : read_ids) {
      m_left_read_ids.insert(read_id);
    }
    m_left.emplace(std::make_pair(a->right_offset, read_ids));

    if (a->matches_reference && !m_options.trace_reference_assemblies) {
      // Discard reference assembly.
      continue;
    }

    calc_score(*a);
    m_output->add(std::move(a));
  }
}

void pair_counter::advance_right() {
  while (!m_right.empty() && m_right.begin()->first < m_cur_offset) {
    auto next_right = m_right.begin();
    assembly_ptr a = std::move(next_right->second);
    m_right.erase(next_right);

    for (const auto& read_id : a->rc_read_ids) {
      auto it = m_right_read_ids.find(read_id);
      CHECK(it != m_right_read_ids.end());
      m_right_read_ids.erase(it);
    }

    a->left_pair_matches =
        find_pair_matches(*a, m_left_read_ids, m_options.forward_pairs_face_inward);
    auto right_offset = a->right_offset;
    m_active.emplace(std::make_pair(right_offset, std::move(a)));
  }
}

void pair_counter::advance_left() {
  while (!m_left.empty() && m_left.begin()->first < (m_cur_offset - m_options.max_pair_distance)) {
    auto next_left = m_left.begin();
    read_id_set read_ids = std::move(next_left->second);
    m_left.erase(next_left);

    for (const auto& read_id : read_ids) {
      auto it = m_left_read_ids.find(read_id);
      CHECK(it != m_left_read_ids.end());
      m_left_read_ids.erase(it);
    }
  }
}

std::vector<uint32_t> pair_counter::find_pair_matches(
    const assembly& a, const std::unordered_multiset<uint32_t, unsalted_hash>& read_ids,
    bool forward) const {
  if (k_pair_counter_debug) {
    std::cout << "Calculating pair matches from " << read_ids.size() << " read ids matching "
              << a.rc_read_ids.size() << " assembly read ids\n";
  }
  std::vector<uint32_t> matches;
  for (const auto& read_id : a.rc_read_ids) {
    if (!m_options.readmap->has_mate(read_id)) {
      continue;
    }
    if (m_options.readmap->get_is_forward(read_id) != forward) {
      continue;
    }
    uint32_t rc_mate_read_id =
        m_options.readmap->get_rev_comp(m_options.readmap->get_mate(read_id));

    if (read_ids.find(rc_mate_read_id) != read_ids.end()) {
      matches.push_back(read_id);
    }
  }
  if (k_pair_counter_debug) {
    std::cout << "Found " << matches.size() << " pair matches\n";
  }
  return matches;
}

}  // namespace variants

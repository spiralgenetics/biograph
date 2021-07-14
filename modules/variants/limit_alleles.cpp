#include "modules/variants/limit_alleles.h"

#include <boost/icl/interval_map.hpp>

namespace variants {

static constexpr int k_dbg = false;

limit_alleles::~limit_alleles() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
  CHECK(m_block_contents.empty());
}

void limit_alleles::on_assembly(assembly_ptr a) {
  aoffset_t left_offset = min(a->left_offset, a->right_offset);
  aoffset_t right_offset = max(a->left_offset, a->right_offset);
  advance_to(left_offset);

  track_left_offset(left_offset);
  m_active.emplace(right_offset, std::move(a));
}

void limit_alleles::advance_to(aoffset_t target) {
  while (m_cur_offset < target) {
    advance_towards(target);
    flush_sorted_to(m_cur_offset);
  }
}

void limit_alleles::advance_towards(aoffset_t target) {
  if (!m_active.empty() && m_active.begin()->first < target) {
    target = m_active.begin()->first;
  }

  m_cur_offset = target;

  while (!m_active.empty() && m_active.begin()->first == m_cur_offset) {
    auto act = m_active.begin();
    m_block_contents.emplace_back(std::move(act->second));
    m_active.erase(act);
  }

  if (m_active.empty() && !m_block_contents.empty()) {
    flush_block_contents();
  }
}

void limit_alleles::flush_block_contents() {
  if (k_dbg) {
    std::cout << "Flushing block at " << m_cur_offset << "\n";
  }
  depths_t depths;
  for (const auto& a : m_block_contents) {
    if (k_dbg) {
      std::cout << "  " << *a << "\n";
    }
    auto i = interval_for_assembly(a);
    depths += std::make_pair(i, size_t(1));
  }
  if (k_dbg) {
    std::cout << "Depths:\n" << depths << "\n";
  }

  if (is_exceeded(depths)) {
    sort_and_limit_block_contents();
  }

  for (auto& a : m_block_contents) {
    untrack_left_offset(min(a->left_offset, a->right_offset));
    sort_and_output(std::move(a));
  }
  m_block_contents.clear();
}

bool limit_alleles::is_exceeded(const depths_t& depths) const {
  for (const auto& i : depths) {
    if (i.second > m_max_alleles) {
      return true;
      break;
    }
  }
  return false;
}

limit_alleles::interval_t limit_alleles::interval_for_assembly(const assembly_ptr& a) {
  // Our interval space includes a position for both interbase and on-base points.
  aoffset_t left_offset = min(a->left_offset, a->right_offset);
  aoffset_t right_offset = max(a->left_offset, a->right_offset);

  if (left_offset == right_offset) {
    // Insert; count for the interbase position between the two bases.
    return interval_t(left_offset * 2, right_offset * 2 + 1);
  } else {
    // Non-insert; don't include interbase positions on either side.
    return interval_t(left_offset * 2 + 1, right_offset * 2);
  }
}

void limit_alleles::sort_and_limit_block_contents() {
  size_t old_size = m_block_contents.size();
  m_block_contents = m_sort_func(std::move(m_block_contents));
  CHECK_EQ(old_size, m_block_contents.size()) << "Some alleles disappeared during sorting?";
  if (k_dbg) {
    std::cout << "Sorting and limiting.\n";
  }

  depths_t depths;
  for (const auto& a : m_block_contents) {
    auto i = interval_for_assembly(a);
    bool limit_this = false;

    if (k_dbg) {
      std::cout << "Checking " << *a << " interval " << i << " against " << depths << "\n";
    }
    auto it = depths.find(i);
    while (it != depths.end() && i.upper() > it->first.lower()) {
      if (it->second >= m_max_alleles) {
        if (k_dbg) {
          std::cout << "Limited because of " << it->first << " with depth " << it->second << "\n";
        }
        limit_this = true;
        break;
      }
      ++it;
    }

    if (limit_this) {
      m_on_limited_func(a);
    } else {
      depths += std::make_pair(i, size_t(1));
      if (k_dbg) {
        std::cout << "New depths: " << depths << "\n";
      }
    }
  }
  CHECK(!is_exceeded(depths));
}

}  // namespace variants

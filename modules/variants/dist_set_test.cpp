#include "modules/variants/dist_set.h"

#include "absl/container/btree_set.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;
using namespace variants;

static constexpr size_t k_mask_bits = dist_set::k_mask_bits;

constexpr uint32_t dist_set_values[] = {0,
                                        1,

                                        k_mask_bits - 1,
                                        k_mask_bits,
                                        k_mask_bits + 1,

                                        2 * k_mask_bits - 1,
                                        2 * k_mask_bits,
                                        2 * k_mask_bits + 1,

                                        3 * k_mask_bits - 1,
                                        3 * k_mask_bits,
                                        3 * k_mask_bits + 1,

                                        4 * k_mask_bits - 1,
                                        4 * k_mask_bits};

static constexpr size_t k_max_dist = 5 * k_mask_bits;

constexpr size_t k_num_dist_set_values = sizeof(dist_set_values) / sizeof(*dist_set_values);
constexpr size_t k_num_small_dist_set_values = 8;

using testing_set = absl::btree_set<int32_t>;

class dist_set_test : public TestWithParam<int /* mask */> {
 public:
  void SetUp() { populate_set(&m_dists, &m_dists_set, GetParam()); }

  void populate_set(dist_set* dists, testing_set* set_dists, unsigned mask) {
    const uint32_t* values = dist_set_values;
    while (mask) {
      if (mask & 1) {
        dists->insert(*values);
        set_dists->insert(*values);
      }
      mask >>= 1;
      ++values;
    }

    ASSERT_THAT(dists_to_set(*dists), ContainerEq(*set_dists));
  }

  testing_set dists_to_set(const dist_set& dists) {
    testing_set result(dists.begin(), dists.end());
    return result;
  }

  testing_set set_add_offset(const testing_set& dists, int to_add) {
    testing_set result;
    for (const auto& dist : dists) {
      result.insert(result.end(), dist + to_add);
    }
    return result;
  }

  testing_set trim_max(testing_set dists, int32_t max_dist, int32_t max_ideal_dist) {
    for (auto it = dists.begin(); it != dists.end(); ++it) {
      if (*it > max_dist) {
        dists.erase(it, dists.end());
        return dists;
      }
      if (*it > max_ideal_dist) {
        ++it;
        dists.erase(it, dists.end());
        return std::move(dists);
      }
    }
    return std::move(dists);
  }

  using closest_table_t = std::map<int32_t, std::pair<int32_t /* first */, int32_t /* last */>>;
  closest_table_t make_closest_table(const dist_set& dists) {
    closest_table_t result;
    for (int target = 0; target < int(k_max_dist); ++target) {
      int32_t best_elem = dists.closest_distance_to(target);

      auto rit = result.find(best_elem);
      if (rit == result.end()) {
        result.emplace(best_elem, std::make_pair(target, target));
      } else {
        CHECK_EQ(rit->second.second, target - 1);
        ++rit->second.second;
      }
    }
    return result;
  }

  closest_table_t set_make_closest_table(const testing_set& dists) {
    closest_table_t result;
    for (int target = 0; target < int(k_max_dist); ++target) {
      int32_t best_elem;
      auto it = dists.upper_bound(target);
      if (it == dists.end()) {
        CHECK(it != dists.begin());
        --it;
        best_elem = *it;
      } else {
        best_elem = *it;
        if (it != dists.begin()) {
          --it;
          if (abs(*it - target) <= abs(best_elem - target)) {
            best_elem = *it;
          }
        }
      }

      auto rit = result.find(best_elem);
      if (rit == result.end()) {
        result.emplace(best_elem, std::make_pair(target, target));
      } else {
        CHECK_EQ(rit->second.second, target - 1);
        ++rit->second.second;
      }
    }
    return result;
  }

 protected:
  dist_set m_dists;
  testing_set m_dists_set;
};

TEST_P(dist_set_test, closest) {
  if (m_dists_set.empty()) {
    return;
  }
  auto actual = make_closest_table(m_dists);
  auto expected = set_make_closest_table(m_dists_set);

  EXPECT_THAT(actual, ContainerEq(expected))
      << "When calculating dist set for: " << PrintToString(m_dists_set)
      << "\nActual: " << PrintToString(actual) << "\nExpected:" << PrintToString(expected);
}

TEST_P(dist_set_test, add_offset) {
  for (size_t offset = 0; offset < k_max_dist; ++offset) {
    auto expected = set_add_offset(m_dists_set, offset);
    auto actual = m_dists.add_offset(offset);
    EXPECT_THAT(dists_to_set(actual), ContainerEq(expected))
        << "\nResult: " << PrintToString(actual) << "\nExpected: " << PrintToString(expected)
        << "\nWhen adding adding offset " << offset << "\nTo: " << PrintToString(m_dists_set);
  }
}

INSTANTIATE_TEST_CASE_P(dist_set_tests, dist_set_test,
                        ::testing::Range(0, 1 << k_num_dist_set_values));

class small_dist_set_test : public dist_set_test {};

TEST_P(small_dist_set_test, add_offset_max) {
  const int32_t max_ideal_dist = std::numeric_limits<int32_t>::max();
  for (size_t offset = 0; offset < k_max_dist; ++offset) {
    for (int32_t max_dist = 0; max_dist < int32_t(k_max_dist); ++max_dist) {
      auto expected = set_add_offset(m_dists_set, offset);
      auto actual = m_dists.add_offset(offset, max_dist, max_ideal_dist);
      auto trimmed_expected = trim_max(expected, max_dist, max_ideal_dist);
      auto trimmed_actual = trim_max(dists_to_set(actual), max_dist, max_ideal_dist);
      if (trimmed_expected != trimmed_actual) {
        EXPECT_THAT(trimmed_actual, ContainerEq(trimmed_expected))
            << "\nResult: " << PrintToString(actual) << "\nExpected: " << PrintToString(expected)
            << "\nWhen adding adding offset " << offset << "\nTo: " << PrintToString(m_dists_set)
            << " with max= " << max_dist << " max ideal=" << max_ideal_dist;
      }
    }
  }
}

TEST_P(small_dist_set_test, add_offset_max_ideal) {
  const int32_t max_dist = std::numeric_limits<int32_t>::max();
  for (size_t offset = 0; offset < k_max_dist; ++offset) {
    for (int32_t max_ideal_dist = 0; max_ideal_dist < int32_t(k_max_dist); ++max_ideal_dist) {
      auto expected = set_add_offset(m_dists_set, offset);
      auto actual = m_dists.add_offset(offset, max_dist, max_ideal_dist);
      auto trimmed_expected = trim_max(expected, max_dist, max_ideal_dist);
      auto trimmed_actual = trim_max(dists_to_set(actual), max_dist, max_ideal_dist);
      if (trimmed_expected != trimmed_actual) {
        EXPECT_THAT(trimmed_actual, ContainerEq(trimmed_expected))
            << "\nResult: " << PrintToString(actual) << "\nExpected: " << PrintToString(expected)
            << "\nWhen adding adding offset " << offset << "\nTo: " << PrintToString(m_dists_set)
            << " with max= " << max_dist << " max ideal=" << max_ideal_dist;
      }
    }
  }
}

INSTANTIATE_TEST_CASE_P(small_dist_set_tests, small_dist_set_test,
                        ::testing::Range(0, 1 << k_num_small_dist_set_values));

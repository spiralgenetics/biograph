#include "modules/variants/read_set.h"

#include "modules/test/set_ops.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;
using namespace variants;

constexpr uint32_t read_id_set_values[] = {0, 1, 1000, 1001, 2000, 2001};
constexpr size_t k_num_read_id_set_values =
    sizeof(read_id_set_values) / sizeof(*read_id_set_values);

class read_id_set_test : public TestWithParam<std::tuple<int /* lhs mask */, int /* rhs mask */>> {
 public:
  void SetUp() {
    populate_set(&m_lhs, &m_lhs_set, std::get<0>(GetParam()));
    populate_set(&m_rhs, &m_rhs_set, std::get<1>(GetParam()));
  }

  void populate_set(read_id_set* reads, std::set<uint32_t>* set_reads, unsigned mask) {
    const uint32_t* values = read_id_set_values;
    while (mask) {
      if (mask & 1) {
        reads->insert(*values);
        set_reads->insert(*values);
      }
      mask >>= 1;
      ++values;
    }

    ASSERT_THAT(read_ids_to_set(*reads), ContainerEq(*set_reads));
  }

  std::set<uint32_t> read_ids_to_set(const read_id_set& reads) {
    std::set<uint32_t> result(reads.begin(), reads.end());
    return result;
  }

 protected:
  read_id_set m_lhs;
  std::set<uint32_t> m_lhs_set;
  read_id_set m_rhs;
  std::set<uint32_t> m_rhs_set;

  read_id_set m_result;
  std::set<uint32_t> m_result_set;
};

TEST_P(read_id_set_test, set_union) {
  m_result = m_lhs | m_rhs;

  std::set_union(m_lhs_set.begin(), m_lhs_set.end(), m_rhs_set.begin(), m_rhs_set.end(),
                 std::inserter(m_result_set, m_result_set.end()));

  EXPECT_THAT(read_ids_to_set(m_result), ContainerEq(m_result_set))
      << "\nLhs: " << PrintToString(m_lhs_set) << "\nRhs: " << PrintToString(m_rhs_set);
}

TEST_P(read_id_set_test, set_difference) {
  m_result = m_lhs - m_rhs;

  std::set_difference(m_lhs_set.begin(), m_lhs_set.end(), m_rhs_set.begin(), m_rhs_set.end(),
                      std::inserter(m_result_set, m_result_set.end()));

  EXPECT_THAT(read_ids_to_set(m_result), ContainerEq(m_result_set))
      << "\nLhs: " << PrintToString(m_lhs_set) << "\nRhs: " << PrintToString(m_rhs_set);
}

TEST_P(read_id_set_test, set_intersection) {
  m_result = m_lhs & m_rhs;

  std::set_intersection(m_lhs_set.begin(), m_lhs_set.end(), m_rhs_set.begin(), m_rhs_set.end(),
                        std::inserter(m_result_set, m_result_set.end()));

  EXPECT_THAT(read_ids_to_set(m_result), ContainerEq(m_result_set))
      << "\nLhs: " << PrintToString(m_lhs_set) << "\nRhs: " << PrintToString(m_rhs_set);
}

INSTANTIATE_TEST_CASE_P(read_id_set_tests, read_id_set_test,
                        ::testing::Combine(  //
                            ::testing::Range(0, 1 << k_num_read_id_set_values),
                            ::testing::Range(0, 1 << k_num_read_id_set_values)));

read_id_set read_id_set_for_elems(const std::vector<uint32_t>& elems) {
  read_id_set new_container;
  new_container.insert(elems.begin(), elems.end());
  return new_container;
}

big_read_id_set big_read_id_set_for_elems(const std::vector<uint32_t>& elems) {
  read_id_set new_container = read_id_set_for_elems(elems);
  big_read_id_set result;
  result |= new_container;
  return result;
}

// Test set operations for read_id_set

struct read_id_set_test_traits
    : set_ops_test_traits_base<read_id_set_test_traits, uint32_t /* read id */, read_id_set> {
  static std::vector<uint32_t> example_elems() {
    return {// First chunk
            0, 1,
            // Second chunk:
            (1 << 10) + 1, (1 << 10) + 2,
            // Third chunk:
            2 << 10, (2 << 10) + 1};
  }

  static read_id_set container_for_elems(const std::vector<uint32_t>& elems) {
    return read_id_set_for_elems(elems);
  }
  static read_id_set rhs_container_for_elems(const std::vector<uint32_t>& elems) {
    return read_id_set_for_elems(elems);
  }
};

INSTANTIATE_TYPED_TEST_CASE_P(ReadIdSet, SetOpsTests, read_id_set_test_traits);

static read_coverage_read_t cov_read_for_read(uint32_t read_id) {
  switch (read_id) {
    case 0:
      return read_coverage_read_t(10, 0, 100);
    case 1:
      return read_coverage_read_t(10, 1, 100);
    case 2:
      return read_coverage_read_t(10, 2, 150);
    case 3:
      return read_coverage_read_t(10, 3, 150);
    case 4:
      return read_coverage_read_t(100, 4, 150);
    case 5:
      return read_coverage_read_t(100, 5, 150);
    default:
      LOG(FATAL) << "Unknown read id: " << read_id;
  }
};

static std::vector<uint32_t> coverage_read_ids(const read_coverage_t& cov) {
  std::vector<uint32_t> result;
  for (const auto& cov_entry : cov.reads()) {
    CHECK(!cov_entry.read_ids.empty());
    for (const auto& read_id : cov_entry.read_ids) {
      result.push_back(read_id);
    }
  }
  return result;
}

read_coverage_t read_coverage_for_elems(const std::vector<uint32_t> read_ids) {
  read_coverage_set cov;

  for (const auto& read_id : read_ids) {
    cov.insert(cov_read_for_read(read_id));
  }

  return cov.build_and_clear(1000 /* assembly len */);
}

// Test set operations for read_coverage_t
struct read_coverage_test_traits
    : set_ops_test_traits_base<read_coverage_test_traits, uint32_t, read_coverage_t> {
  static std::vector<uint32_t> example_elems() { return {0, 1, 2, 3}; }

  static read_coverage_t container_for_elems(const std::vector<uint32_t>& elems) {
    return read_coverage_for_elems(elems);
  }
  static read_coverage_t rhs_container_for_elems(const std::vector<uint32_t>& elems) {
    return read_coverage_for_elems(elems);
  }

  static std::vector<uint32_t> elems_in_container(const read_coverage_t& cov) {
    return coverage_read_ids(cov);
  }
  static std::vector<uint32_t> rhs_elems_in_container(const read_coverage_t& cov) {
    return coverage_read_ids(cov);
  }
};

INSTANTIATE_TYPED_TEST_CASE_P(ReadCoverage, SetOpsTests, read_coverage_test_traits);

// Test set operations for read_coverage_t versus read_id_set
struct read_coverage_vs_ids_test_traits
    : set_ops_test_traits_base<read_coverage_vs_ids_test_traits, uint32_t, read_coverage_t,
                               read_id_set> {
  static std::vector<uint32_t> example_elems() { return {0, 1, 2, 3, 4, 5}; }

  static read_coverage_t container_for_elems(const std::vector<uint32_t>& elems) {
    return read_coverage_for_elems(elems);
  }
  static read_id_set rhs_container_for_elems(const std::vector<uint32_t>& elems) {
    return read_id_set_test_traits::container_for_elems(elems);
  }

  static std::vector<uint32_t> elems_in_container(const read_coverage_t& cov) {
    return coverage_read_ids(cov);
  }
  static std::vector<uint32_t> rhs_elems_in_container(const read_id_set& ids) {
    return read_id_set_test_traits::elems_in_container(ids);
  }
};

namespace variants {

// read_coverage_set | read_id_set is nonsensical since we wouldn't
// know where to align the added reads, so fake it for the test.
read_coverage_t operator|(const read_coverage_t& lhs, const read_id_set& rhs_ids) {
  read_coverage_t rhs = read_coverage_vs_ids_test_traits::container_for_elems(
      read_coverage_vs_ids_test_traits::rhs_elems_in_container(rhs_ids));
  return lhs | rhs;
}
read_coverage_t& operator|=(read_coverage_t& lhs, const read_id_set& rhs_ids) {
  read_coverage_t rhs = read_coverage_vs_ids_test_traits::container_for_elems(
      read_coverage_vs_ids_test_traits::rhs_elems_in_container(rhs_ids));
  return lhs |= rhs;
}

}  // namespace variants

INSTANTIATE_TYPED_TEST_CASE_P(ReadCoverageVsIds, SetOpsTests, read_coverage_vs_ids_test_traits);

// Test set operations for big_read_id_set versus read_id_set
struct big_ids_vs_ids_test_traits
    : set_ops_test_traits_base<read_coverage_vs_ids_test_traits, uint32_t, big_read_id_set,
                               read_id_set> {
  static std::vector<uint32_t> example_elems() { return {0, 1, 2, 3, 4, 5}; }

  static big_read_id_set container_for_elems(const std::vector<uint32_t>& elems) {
    return big_read_id_set_for_elems(elems);
  }
  static read_id_set rhs_container_for_elems(const std::vector<uint32_t>& elems) {
    return read_id_set_for_elems(elems);
  }
};

INSTANTIATE_TYPED_TEST_CASE_P(BigSetVsSet, SetOpsTests, big_ids_vs_ids_test_traits);

// Test set operations for big_read_id_set versus read_id_set
struct ids_vs_big_ids_test_traits
    : set_ops_test_traits_base<read_coverage_vs_ids_test_traits, uint32_t, read_id_set,
                               big_read_id_set> {
  static std::vector<uint32_t> example_elems() { return {0, 1, 2, 3, 4, 5}; }

  static read_id_set container_for_elems(const std::vector<uint32_t>& elems) {
    return read_id_set_for_elems(elems);
  }
  static big_read_id_set rhs_container_for_elems(const std::vector<uint32_t>& elems) {
    return big_read_id_set_for_elems(elems);
  }
};

INSTANTIATE_TYPED_TEST_CASE_P(SetVsBigSet, SetOpsTests, ids_vs_big_ids_test_traits);

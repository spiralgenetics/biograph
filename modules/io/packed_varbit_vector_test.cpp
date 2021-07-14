#include "modules/io/packed_varbit_vector.h"
#include "modules/io/log.h"
#include "modules/io/spiral_file_mem.h"

#include <random>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;

TEST(packed_varbit_vector_test, sizes) {
  EXPECT_EQ(0, packed_varbit_vector::bits_for_value(0));

  EXPECT_EQ(0, packed_varbit_vector::calc_size(0 /* elems */, 0 /* max value */));
  EXPECT_EQ(0, packed_varbit_vector::calc_size(1 /* elems */, 0 /* max value */));
  EXPECT_EQ(0, packed_varbit_vector::calc_size(1000 /* elems */, 0 /* max value */));

  EXPECT_EQ(1, packed_varbit_vector::bits_for_value(1));

  EXPECT_EQ(8, packed_varbit_vector::calc_size(1 /* elems */, 1 /* max value */));
  EXPECT_EQ(8, packed_varbit_vector::calc_size(1 /* elems */, 2 /* max value */));
  EXPECT_EQ(8, packed_varbit_vector::calc_size(1 /* elems */, 3 /* max value */));

  EXPECT_EQ(8, packed_varbit_vector::calc_size(63 /* elems */, 1 /* max value */));
  EXPECT_EQ(8, packed_varbit_vector::calc_size(64 /* elems */, 1 /* max value */));
  EXPECT_EQ(16, packed_varbit_vector::calc_size(65 /* elems */, 1 /* max value */));

  EXPECT_EQ(2, packed_varbit_vector::bits_for_value(2));

  EXPECT_EQ(8, packed_varbit_vector::calc_size(32 /* elems */, 2 /* max value */));
  EXPECT_EQ(16, packed_varbit_vector::calc_size(33 /* elems */, 2 /* max value */));

  EXPECT_EQ(2, packed_varbit_vector::bits_for_value(3));
  EXPECT_EQ(3, packed_varbit_vector::bits_for_value(4));

  EXPECT_EQ(8, packed_varbit_vector::calc_size(21 /* elems */, 4 /* max value */));
  EXPECT_EQ(16, packed_varbit_vector::calc_size(22 /* elems */, 4 /* max value */));

  uint64_t max63 = std::numeric_limits<uint64_t>::max() >> 1;
  EXPECT_EQ(8, packed_varbit_vector::calc_size(1 /* elems */, max63));
  EXPECT_EQ(16, packed_varbit_vector::calc_size(2 /* elems */, max63));
  EXPECT_EQ(8 * 63, packed_varbit_vector::calc_size(63 /* elems */, max63));
  EXPECT_EQ(8 * 63, packed_varbit_vector::calc_size(64 /* elems */, max63));
  EXPECT_EQ(8 * 64, packed_varbit_vector::calc_size(65 /* elems */, max63));

  uint64_t max64 = std::numeric_limits<uint64_t>::max();
  EXPECT_EQ(8, packed_varbit_vector::calc_size(1 /* elems */, max64));
  EXPECT_EQ(16, packed_varbit_vector::calc_size(2 /* elems */, max64));
}

class packed_varbit_vector_test_p
    : public TestWithParam<std::pair<size_t /* number of elems */, size_t /* max value */>> {
 public:
  packed_varbit_vector_test_p() { std::tie(m_num_elems, m_max_value) = GetParam(); }

 protected:
  size_t m_num_elems;
  size_t m_max_value;
};

TEST_P(packed_varbit_vector_test_p, basic) {
  std::random_device rand_dev;
  size_t seed = rand_dev();
  SPLOG("Using seed %ld", seed);
  std::mt19937_64 rand_source1(seed);
  std::mt19937_64 rand_source2(seed);
  std::uniform_int_distribution<uint64_t> rand_value(0, m_max_value);

  spiral_file_create_mem c;
  {
    mutable_packed_varbit_vector vec(c.create(), m_num_elems, m_max_value);

    for (size_t i = 0; i < m_num_elems; ++i) {
      vec.set(i, rand_value(rand_source1));
    }
    for (size_t i = 0; i < m_num_elems; ++i) {
      EXPECT_EQ(rand_value(rand_source2), vec.get(i)) << "i: " << i;
    }
    // Do it again, to make sure updating works.
    for (size_t i = 0; i < m_num_elems; ++i) {
      vec.set(i, rand_value(rand_source1));
    }
  }

  spiral_file_open_mem o(c.close());
  packed_varbit_vector vec(o.open());
  EXPECT_EQ(m_num_elems, vec.size());
  EXPECT_EQ(m_max_value, vec.max_value());
  for (size_t i = 0; i < m_num_elems; ++i) {
    EXPECT_EQ(rand_value(rand_source2), vec.get(i)) << "i: " << i;
  }
}

INSTANTIATE_TEST_CASE_P(
    packed_varbit_vector_tests, packed_varbit_vector_test_p,
    ::testing::ValuesIn(std::vector<std::pair<size_t /* elem count */, size_t /* max value */>>{
        {0, 0},
        {1, 0},
        {100, 0},
        {1, 1},
        {63, 1},
        {64, 1},
        {65, 1},
        {63, 2},
        {64, 2},
        {65, 2},
        {21, 4},
        {21, 4},
        {22, 4},
        {22, 5},
        // 63 bits per element
        {63, 1ULL << 62},
        {64, 1ULL << 62},
        {65, 1ULL << 62},
        {66, 1ULL << 62},
        // 64 bits per element
        {63, 3ULL << 62},
        {64, 3ULL << 62},
        {65, 3ULL << 62},
        {66, 3ULL << 62},
    }));

class packed_varbit_vector_thread_test_p
    : public TestWithParam<std::pair<size_t /* number of elems */, size_t /* max value */>> {
 public:
  packed_varbit_vector_thread_test_p() { std::tie(m_num_elems, m_max_value) = GetParam(); }

 protected:
  size_t m_num_elems;
  size_t m_max_value;
};

TEST_P(packed_varbit_vector_thread_test_p, threads) {
  const size_t k_num_threads = 8;
  const size_t k_iter_count = 100000;

  spiral_file_create_mem c;
  {
    mutable_packed_varbit_vector pvec(c.create(), m_num_elems, m_max_value);

    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < k_num_threads; i++) {
      auto run_thread = [this, &pvec](size_t thread_id) {
        std::vector<size_t> expected;
        std::random_device rand_dev;
        size_t seed = rand_dev();
        std::mt19937_64 rand_source(seed);
        size_t num_this_thread = m_num_elems / k_num_threads;
        if (k_num_threads * num_this_thread + thread_id < m_num_elems) {
          ++num_this_thread;
        }
        expected.resize(m_num_elems / k_num_threads, 0);
        std::uniform_int_distribution<size_t> elem_idx(0, expected.size() - 1);
        std::uniform_int_distribution<uint64_t> rand_value(0, m_max_value);

        for (size_t count = 0; count < k_iter_count; ++count) {
          size_t idx = elem_idx(rand_source);

          size_t new_val = rand_value(rand_source);
          size_t old_val = pvec.get(idx * k_num_threads + thread_id);
          pvec.set(idx * k_num_threads + thread_id, new_val);
          ASSERT_EQ(expected[idx], old_val);
          expected[idx] = new_val;
        }
      };
      futures.emplace_back(std::async(std::launch::async, run_thread, i));
    }

    for (std::future<void>& future : futures) {
      future.get();
    }
  }
}

INSTANTIATE_TEST_CASE_P(
    varbit_vector_thread_tests, packed_varbit_vector_thread_test_p,
    ::testing::ValuesIn(std::vector<std::pair<size_t /* elem count */, size_t /* max value */>>{
        {1000, 0},
        {1027, 1},
        // 5 bits per elem
        {1029, (1ULL << 5) - 1},
        {1030, (1ULL << 5) - 1},
        {1031, (1ULL << 5) - 1},
        // 63 bits per element
        {1007, (~uint64_t(0)) >> 1},
        {1008, (~uint64_t(0)) >> 1},
        {1009, (~uint64_t(0)) >> 1},
        // 64 bits per element
        {1047, (~uint64_t(0))}}));

class packed_varbit_vector_uint8_test_p : public TestWithParam<int /* size of membuf */> {
 public:
  packed_varbit_vector_uint8_test_p() { m_num_elems = GetParam(); }

 protected:
  size_t m_num_elems;
};

INSTANTIATE_TEST_CASE_P(varbit_vector_membuf_tests, packed_varbit_vector_uint8_test_p,
                        ::testing::Values(0,   1,   2,   3,   4,   5,   6,   7,   8,   9,
                                             900, 901, 902, 903, 904, 905, 906, 907, 908, 909));

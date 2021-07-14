#include "modules/io/sparse_multi.h"
#include "modules/io/spiral_file_mem.h"

#include <random>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;

namespace {

struct test_case_part {
  size_t src_count = 0;
  size_t dest_count = 0;
  std::string desc;
  test_case_part(size_t src, size_t dest, const std::string& d)
      : src_count(src), dest_count(dest), desc(d){};
};

std::ostream& operator<<(std::ostream& os, const test_case_part& part) {
  return os << part.desc;
}

test_case_part multi(size_t dest_count) {
  return test_case_part{1, dest_count,
                        "multi(" + std::to_string(dest_count) + ")"};
}

test_case_part skip(size_t src_count = 1) {
  return test_case_part{src_count, 0,
                        "skip(" + std::to_string(src_count) + ")"};
}

test_case_part single() { return test_case_part{1, 1, "single"}; }

test_case_part singles(size_t src_count) {
  return test_case_part{src_count, 1,
                        "singles(" + std::to_string(src_count) + ")"};
}

std::vector<std::vector<test_case_part>> test_cases{
    {},
    {skip()},
    {skip((2 << 16) - 1)},
    {skip(2 << 16)},
    {skip((2 << 16) + 1)},

    {single()},
    {single(), skip()},
    {skip(), single()},
    {multi(2)},
    {multi(2), skip()},
    {skip(), multi(2)},

    {skip((1 << 16) - 1), multi(2)},
    {skip(1 << 16), multi(2)},
    {skip((1 << 16) + 1), multi(2)},

    {multi(2), skip((1 << 16) - 1)},
    {multi(2), skip(1 << 16)},
    {multi(2), skip((1 << 16) + 1)},

    {skip((2 << 16) - 1), multi(2)},
    {skip(2 << 16), multi(2)},
    {skip((2 << 16) + 1), multi(2)},

    {multi(2), skip((2 << 16) - 1)},
    {multi(2), skip(2 << 16)},
    {multi(2), skip((2 << 16) + 1)},

    {skip((2 << 16) - 1), multi(2), skip((2 << 16) - 1), multi(2),
     skip((2 << 16) - 1), multi(2)},
    {skip(2 << 16), multi(2), skip(2 << 16), multi(2), skip(2 << 16), multi(2)},
    {skip((2 << 16) + 1), multi(2), skip((2 << 16) + 1), multi(2),
     skip((2 << 16) + 1), multi(2)},

    {singles((1 << 16) - 1)},
    {singles(1 << 16)},
    {singles((1 << 16) + 1)},
};

}  // namespace

class sparse_multi_test
    : public testing::TestWithParam<std::vector<test_case_part>> {};

TEST_P(sparse_multi_test, verify_range) {
  const std::vector<test_case_part>& parts = GetParam();

  std::cerr << ::testing::PrintToString(parts) << "\n";

  size_t num_srcs = 0, num_dsts = 0;

  for (const test_case_part& part : parts) {
    num_srcs += part.src_count;
    num_dsts += part.dest_count * part.src_count;
  }

  spiral_file_create_mem c;
  sparse_multi_builder b(c.create(), num_srcs, num_dsts);

  size_t cur_src = 0;
  size_t cur_dst = 0;
  for (const test_case_part& part : parts) {
    for (size_t i = 0; i < part.src_count; i++) {
      for (size_t j = 0; j < part.dest_count; j++) {
        EXPECT_EQ(cur_dst, b.add(cur_src));
        cur_dst++;
      }
      cur_src++;
    }
  }

  EXPECT_EQ(cur_src, num_srcs);
  EXPECT_EQ(cur_dst, num_dsts);

  auto multi = b.finalize();

  std::random_device rand_dev;
  std::mt19937 rand_source(rand_dev());
  std::uniform_int_distribution<int64_t> rand_src(0, num_srcs);
  for (size_t count = 0; count < 1000; ++count) {
    size_t src1 = rand_src(rand_source);
    size_t src2 = rand_src(rand_source);
    if (src1 > src2) {
      std::swap(src1, src2);
    }

    boost::optional<size_t> expected_start, expected_limit;
    for (size_t src = src1; src != src2; src++) {
      auto range = multi->lookup(src);
      if (range.first == range.second) {
        continue;
      }
      if (!expected_start) {
        expected_start.emplace(range.first);
      }
      if (expected_limit) {
        EXPECT_EQ(range.first, *expected_limit);
      }
      expected_limit.emplace(range.second);
    }

    auto range = multi->lookup_range(src1, src2);
    if (!expected_start) {
      EXPECT_FALSE(expected_limit) << " range: " << src1 << " to " << src2
                                   << " actual: " << range.first
                                   << " to: " << range.second;
      EXPECT_EQ(range.first, range.second) << " range: " << src1 << " to "
                                           << src2 << " actual: " << range.first
                                           << " to: " << range.second;
    } else {
      ASSERT_TRUE(expected_limit) << " range: " << src1 << " to " << src2
                                  << " actual: " << range.first
                                  << " to: " << range.second;
      EXPECT_EQ(*expected_start, range.first)
          << " range: " << src1 << " to " << src2 << " actual: " << range.first
          << " to: " << range.second << " expected: " << *expected_start
          << " to: " << *expected_limit;
      EXPECT_EQ(*expected_limit, range.second)
          << " range: " << src1 << " to " << src2 << " actual: " << range.first
          << " to: " << range.second << " expected: " << *expected_start
          << " to: " << *expected_limit;
    }
  }
}

TEST_P(sparse_multi_test, verify) {
  const std::vector<test_case_part>& parts = GetParam();

  std::cerr << ::testing::PrintToString(parts) << "\n";

  size_t num_srcs = 0, num_dsts = 0;

  for (const test_case_part& part : parts) {
    num_srcs += part.src_count;
    num_dsts += part.dest_count * part.src_count;
  }

  spiral_file_create_mem c;
  sparse_multi_builder b(c.create(), num_srcs, num_dsts);

  size_t cur_src = 0;
  size_t cur_dst = 0;
  for (const test_case_part& part : parts) {
    for (size_t i = 0; i < part.src_count; i++) {
      for (size_t j = 0; j < part.dest_count; j++) {
        EXPECT_EQ(cur_dst, b.add(cur_src));
        cur_dst++;
      }
      cur_src++;
    }
  }

  EXPECT_EQ(cur_src, num_srcs);
  EXPECT_EQ(cur_dst, num_dsts);

  auto multi = b.finalize();

  cur_src = cur_dst = 0;

  for (const test_case_part& part : parts) {
    for (size_t i = 0; i < part.src_count; i++) {
      std::vector<size_t> expected_range;
      for (size_t j = 0; j < part.dest_count; j++) {
        expected_range.push_back(cur_dst);
        EXPECT_EQ(cur_src, multi->reverse_lookup(cur_dst));
        cur_dst++;
      }
      auto range = multi->lookup(cur_src);
      std::vector<size_t> actual_range;
      for (size_t j = range.first; j != range.second; j++) {
        actual_range.push_back(j);
      }
      EXPECT_EQ(expected_range, actual_range);
      cur_src++;
    }
  }

  // Next, verify that the deserialized version is OK.
  spiral_file_mem_storage encoded = c.close();
  spiral_file_open_mem o(encoded);
  sparse_multi decoded(o.open());

  cur_src = cur_dst = 0;

  for (const test_case_part& part : parts) {
    for (size_t i = 0; i < part.src_count; i++) {
      EXPECT_EQ(cur_dst, decoded.lookup_lower_bound(cur_src));
      std::vector<size_t> expected_dests;
      for (size_t j = 0; j < part.dest_count; j++) {
        expected_dests.push_back(cur_dst);
        EXPECT_EQ(cur_src, decoded.reverse_lookup(cur_dst));
        if (j == 0) {
          EXPECT_TRUE(decoded.dest_is_first_in_group(cur_dst));
        } else if (decoded.dest_is_first_in_group(cur_dst)) {
          EXPECT_FALSE(decoded.dest_is_first_in_group(cur_dst));
        }
        cur_dst++;
      }
      auto range = decoded.lookup(cur_src);
      std::vector<size_t> actual_dests;
      for (size_t j = range.first; j != range.second; j++) {
        actual_dests.push_back(j);
      }
      EXPECT_EQ(expected_dests, actual_dests);

      cur_src++;
    }
  }

  cur_src = cur_dst = 0;
  std::vector<std::pair<size_t, std::pair<size_t, size_t>>> expected_ranges;
  sparse_multi::iterator cur_it = decoded.begin();
  for (const test_case_part& part : parts) {
    for (size_t i = 0; i < part.src_count; i++) {
      if (i == 0 || i + 1 == part.src_count) {
        EXPECT_EQ(cur_it, decoded.iterator_at_source(cur_src));
      }
      if (part.dest_count) {
        EXPECT_EQ(cur_src, (*cur_it).first);
        EXPECT_EQ(cur_dst, (*cur_it).second.first);
        size_t dest_end = cur_dst + part.dest_count;
        EXPECT_EQ(dest_end, (*cur_it).second.second);
        expected_ranges.emplace_back(
            std::make_pair(cur_src, std::make_pair(cur_dst, dest_end)));
        cur_dst = dest_end;
        ++cur_it;
      }
      cur_src++;
    }
  }
  EXPECT_EQ(decoded.end(), cur_it);

  // And that using the range for loop works.
  auto expected_range_it = expected_ranges.begin();
  for (const auto& actual : decoded) {
    ASSERT_TRUE(expected_range_it != expected_ranges.end());
    EXPECT_EQ(*expected_range_it, actual);
    expected_range_it++;
  }
}

// Build using old "readmap" gross/fine buffers.
TEST_P(sparse_multi_test, build_from_old_format) {
  const std::vector<test_case_part>& parts = GetParam();

  std::cerr << ::testing::PrintToString(parts) << "\n";

  size_t num_srcs = 0, num_dsts = 0;

  for (const test_case_part& part : parts) {
    num_srcs += part.src_count;
    num_dsts += part.dest_count * part.src_count;
  }

  if (num_srcs == 0) {
    // Don't test size 0 seqsets for old format. */
    return;
  }

  // Calculate gross and fine tables.
  size_t last_elem_hi = (num_srcs - 1) >> 16;
  size_t n_gross_ids = last_elem_hi + 2;
  std::vector<uint32_t> gross(n_gross_ids);
  std::vector<uint16_t> fine(num_dsts);

  uint32_t* gross_ptr = gross.data();
  uint16_t* fine_ptr = fine.data();

  size_t cur_src = 0, cur_dst = 0;
  size_t last_source_hi = 0;
  for (const test_case_part& part : parts) {
    for (size_t i = 0; i < part.src_count; i++) {
      uint32_t source_hi = cur_src >> 16;
      for (size_t j = 0; j < part.dest_count; j++) {
        if (gross_ptr == gross.data()) {
          *gross_ptr++ = cur_dst;
        }
        while (source_hi != last_source_hi) {
          *gross_ptr++ = cur_dst;
          last_source_hi++;
        }

        *fine_ptr++ = cur_src & 0xFFFF;
        cur_dst++;
      }
      cur_src++;
    }
  }

  last_elem_hi++;
  while (last_elem_hi != last_source_hi) {
    CHECK_LT(last_source_hi, last_elem_hi);
    *gross_ptr++ = num_dsts;
    last_source_hi++;
  }

  ASSERT_EQ(cur_src, num_srcs);
  ASSERT_EQ(cur_dst, num_dsts);

  spiral_file_create_mem c;
  sparse_multi_builder b(c.create(), num_srcs, num_dsts);
  b.build_from_old_format(reinterpret_cast<char*>(gross.data()),
                          reinterpret_cast<char*>(fine.data()));

  auto multi = b.finalize();

  cur_src = cur_dst = 0;
  for (const test_case_part& part : parts) {
    for (size_t i = 0; i < part.src_count; i++) {
      std::vector<size_t> expected_range;
      for (size_t j = 0; j < part.dest_count; j++) {
        expected_range.push_back(cur_dst);
        EXPECT_EQ(cur_src, multi->reverse_lookup(cur_dst));
        cur_dst++;
      }
      auto range = multi->lookup(cur_src);
      std::vector<size_t> actual_range;
      for (size_t j = range.first; j != range.second; j++) {
        actual_range.push_back(j);
      }
      EXPECT_EQ(expected_range, actual_range);
      cur_src++;
    }
  }
}

INSTANTIATE_TEST_CASE_P(sparse_multi_tests, sparse_multi_test,
                        ValuesIn(test_cases));

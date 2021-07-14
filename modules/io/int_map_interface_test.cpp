#include "modules/io/int_map_interface.h"
#include "modules/io/log.h"
#include "modules/io/packed_vector.h"

#include "gtest/gtest.h"

TEST(less_than_search_test, search_backward) {
  mutable_packed_vector<size_t, 32> pvec(64 * 64 + 37, "less_than_search_test:search_backward");
  static constexpr size_t k_offset = 123;
  for (size_t i = 0; i != pvec.size(); ++i) {
    pvec[i] = i + k_offset - 1;
  }
  less_than_search search(&pvec);

  // Interesting locations to test
  std::set<int> int_locs;
  int size = pvec.size();
  int k_factor1 = less_than_search::get_factor1();
  int k_factor2 = less_than_search::get_factor2();
  for (int j : {1 - k_factor1, -1, 0, 1, k_factor1 - 1}) {
    for (int i : {0, k_factor1 * k_factor2, 2 * k_factor1 * k_factor2, size,
                  (size / k_factor1) * (k_factor1),
                  (size / (k_factor1 * k_factor2)) * k_factor1 * k_factor2}) {
      int loc = i + j;
      if (loc < 0 || loc >= int(pvec.size())) {
        continue;
      }
      int_locs.insert(loc);
    }
  }

  for (size_t loc : int_locs) {
    for (size_t i = 0; i != pvec.size(); ++i) {
      if (loc > i) {
        EXPECT_EQ(search.next_backward_lt(loc, i + k_offset), i) << " loc: " << loc << " i: " << i;
        EXPECT_EQ(search.next_forward_lt(loc, i + k_offset), size)
            << " loc: " << loc << " i: " << i;
      } else {
        EXPECT_EQ(search.next_backward_lt(loc, i + k_offset), loc)
            << " loc: " << loc << " i: " << i;
        EXPECT_EQ(search.next_forward_lt(loc, i + k_offset), loc) << " loc: " << loc << " i: " << i;
      }
    }
  }
}

TEST(less_than_search_test, search_forward) {
  mutable_packed_vector<size_t, 32> pvec(64 * 64 + 37, "less_than_search_test:search_backward");
  static constexpr size_t k_offset = 64 * 64 + 37 + 123;
  for (size_t i = 0; i != pvec.size(); ++i) {
    pvec[i] = k_offset - i - 1;
  }
  less_than_search search(&pvec);

  // Interesting locations to test
  std::set<int> int_locs;
  int size = pvec.size();
  int k_factor1 = less_than_search::get_factor1();
  int k_factor2 = less_than_search::get_factor2();
  for (int j : {1 - k_factor1, -1, 0, 1, k_factor1 - 1}) {
    for (int i : {0, k_factor1 * k_factor2, 2 * k_factor1 * k_factor2, size,
                  (size / k_factor1) * (k_factor1),
                  (size / (k_factor1 * k_factor2)) * k_factor1 * k_factor2}) {
      int loc = i + j;
      if (loc < 0 || loc >= int(pvec.size())) {
        continue;
      }
      int_locs.insert(loc);
    }
  }

  for (size_t loc : int_locs) {
    for (size_t i = 0; i != pvec.size(); ++i) {
      if (loc < i) {
        EXPECT_EQ(search.next_forward_lt(loc, k_offset - i), i) << " loc: " << loc << " i: " << i;
        EXPECT_EQ(search.next_backward_lt(loc, k_offset - i), 0) << " loc: " << loc << " i: " << i;
      } else {
        EXPECT_EQ(search.next_forward_lt(loc, k_offset - i), loc) << " loc: " << loc << " i: " << i;
        EXPECT_EQ(search.next_backward_lt(loc, k_offset - i), loc) << " loc: " << loc << " i: " << i;
      }
    }
  }
}

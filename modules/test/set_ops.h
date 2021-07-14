#pragma once

#include "base/base.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <vector>

#include <boost/optional.hpp>

// Test suite to test set operations on an arbitrary type to make sure they're consistent.
//
// This tests the following operations:
//
//   lhs | rhs:  set union
//   lhs & rhs:  set intersection
//   lhs - rhs:  set difference
//
//   lhs |=,&=,-= rhs: Union, intersection and difference in place.
//
//   empty():  Returns true if empty.
//
// This is parameterized to a "Traits" class, which must supply
// information on the container class being tested.
//
// Here is an example minimal "Traits" class for a "MySet" class which contains "int"s.
//
// struct my_set_ops_test_traits : set_ops_test_traits_base<my_set_ops_test_traits,
//      int /* element type */, my_set /* container type */> {
//   // Function which returns some example elements of type elem_t.
//   static std::vector<int> example_elems() {
//     return {1, 2, 3, 4};
//   };
// };
//
// Note that there are additional traits available with defaults in
// my_set_ops_test_traits_base that can be overriden.
//
// Set operations can also be tested between two different conainer types.  For example,
// to test "lhs | rhs", "lhs & rhs", etc:
//
// struct my_lhs_rhs_set_ops_test_traits : set_ops_test_traits_base<my_set_ops_test_traits,
//      int /* element type */,
//      lhs_set /* container type on left hand side of operators */,
//      rhs_set /* container type o  right hand side of operators */> {
//   // Function which returns some example elements of type elem_t.
//   static std::vector<int> example_elems() {
//     return {1, 2, 3, 4};
//   };
// };
//
template <typename Traits, typename ElementType, typename ContainerType,
          typename RhsContainerType = ContainerType>
struct set_ops_test_traits_base {
  using container_t = ContainerType;
  using rhs_container_t = RhsContainerType;
  using elem_t = ElementType;

  // Gather elements present in the container into a vector.
  static std::vector<elem_t> elems_in_container(const container_t& container) {
    std::vector<elem_t> elems;
    for (const auto& elem : container) {
      elems.push_back(elem);
    }
    return elems;
  }

  // Gather elements present in the container into a vector.
  static std::vector<elem_t> rhs_elems_in_container(const rhs_container_t& container) {
    std::vector<elem_t> elems;
    for (const auto& elem : container) {
      elems.push_back(elem);
    }
    return elems;
  }

  // Create a container given a vector of elements.
  static container_t container_for_elems(const std::vector<elem_t>& elems) {
    container_t new_container(elems.begin(), elems.end());
    return new_container;
  }
  static rhs_container_t rhs_container_for_elems(const std::vector<elem_t>& elems) {
    rhs_container_t new_container(elems.begin(), elems.end());
    return new_container;
  }

  // Returns true if the element lhs should come before the element rhs in the set.
  static bool elem_less_than(const elem_t& lhs, const elem_t& rhs) { return lhs < rhs; }

  // Returns a sampling of element sets that should be tested against each other
  static std::vector<std::vector<elem_t>> example_element_sets() {
    std::vector<std::vector<elem_t>> examples;

    std::vector<elem_t> all_elems = Traits::example_elems();
    CHECK_LE(all_elems.size(), 10) << "Too many example elements to generate all combinations";

    for (size_t i = 0; i < 1ULL << all_elems.size(); ++i) {
      std::vector<elem_t> this_example;
      size_t include_left = i;
      for (const auto& elem : all_elems) {
        if (include_left & 1) {
          this_example.push_back(elem);
        }
        include_left >>= 1;
      }
      CHECK_EQ(include_left, 0);

      examples.emplace_back(std::move(this_example));
    }
    return examples;
  }
};

template <typename Traits>
class SetOpsTests : public ::testing::Test {
  using container_t = typename Traits::container_t;
  using elem_t = typename Traits::elem_t;
  using elem_vec_t = std::vector<elem_t>;

 protected:
  // Verifies container is consistent.
  void check_container(const container_t& container) {
    boost::optional<elem_t> last_elem;
    for (const auto& elem : Traits::elems_in_container(container)) {
      if (last_elem) {
        EXPECT_TRUE(Traits::elem_less_than(*last_elem, elem));
      }
      last_elem.emplace(elem);
    }

    if (last_elem) {
      EXPECT_FALSE(container.empty());
    } else {
      EXPECT_TRUE(container.empty());
    }
  }

  elem_vec_t set_union(const elem_vec_t& lhs, const elem_vec_t& rhs) {
    elem_vec_t result;
    std::set_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::back_inserter(result),
                   Traits::elem_less_than);
    return result;
  }
  elem_vec_t set_difference(const elem_vec_t& lhs, const elem_vec_t& rhs) {
    elem_vec_t result;
    std::set_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::back_inserter(result),
                        Traits::elem_less_than);
    return result;
  }
  elem_vec_t set_intersection(const elem_vec_t& lhs, const elem_vec_t& rhs) {
    elem_vec_t result;
    std::set_intersection(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                          std::back_inserter(result), Traits::elem_less_than);
    return result;
  }

  void check_create_and_access(const elem_vec_t& orig_elems, const elem_vec_t& expected_elems) {
    container_t container = Traits::container_for_elems(orig_elems);
    check_container(container);
    EXPECT_EQ(expected_elems.empty(), container.empty());
    elem_vec_t actual_elems = Traits::elems_in_container(container);

    EXPECT_THAT(actual_elems, testing::ContainerEq(expected_elems));
  }

  void for_all_example_pairs(
      const std::function<void(const elem_vec_t& lhs, const elem_vec_t& rhs)>& check_pair_f) {
    std::vector<elem_vec_t> examples = Traits::example_element_sets();
    for (const elem_vec_t& lhs : examples) {
      for (const elem_vec_t& rhs : examples) {
        check_pair_f(lhs, rhs);
      }
    }
  }
};

TYPED_TEST_CASE_P(SetOpsTests);

TYPED_TEST_P(SetOpsTests, Empty) {
  using container_t = typename TypeParam::container_t;

  container_t empty_container;
  this->check_container(empty_container);
  EXPECT_TRUE(empty_container.empty());

  container_t empty_container2 = TypeParam::container_for_elems({});
  this->check_container(empty_container2);
  EXPECT_TRUE(empty_container2.empty());
}

TYPED_TEST_P(SetOpsTests, CreateAndAccess) {
  using elem_vec_t = std::vector<typename TypeParam::elem_t>;

  elem_vec_t example_elems = TypeParam::example_elems();
  this->check_create_and_access(example_elems, example_elems);

  elem_vec_t rev_elems = example_elems;
  std::reverse(rev_elems.begin(), rev_elems.end());
  this->check_create_and_access(rev_elems, example_elems);

  elem_vec_t dup_elems = example_elems;
  dup_elems.insert(dup_elems.end(), example_elems.begin(), example_elems.end());
  this->check_create_and_access(dup_elems, example_elems);

  while (!example_elems.empty()) {
    example_elems.pop_back();
    this->check_create_and_access(example_elems, example_elems);
  }
}

TYPED_TEST_P(SetOpsTests, SetUnion) {
  using elem_vec_t = std::vector<typename TypeParam::elem_t>;
  using container_t = typename TypeParam::container_t;
  using rhs_container_t = typename TypeParam::rhs_container_t;

  this->for_all_example_pairs([this](const elem_vec_t& lhs, const elem_vec_t& rhs) {
    container_t lhs_c = TypeParam::container_for_elems(lhs);
    rhs_container_t rhs_c = TypeParam::rhs_container_for_elems(rhs);
    elem_vec_t expected = this->set_union(lhs, rhs);

    elem_vec_t actual = TypeParam::elems_in_container(lhs_c | rhs_c);
    EXPECT_THAT(actual, testing::ContainerEq(expected))
        << "From " << testing::PrintToString(lhs) << " | " << testing::PrintToString(rhs);

    auto* in_place_result = &(lhs_c |= rhs_c);
    EXPECT_EQ(in_place_result, &lhs_c);
    elem_vec_t actual_in_place = TypeParam::elems_in_container(lhs_c);
    EXPECT_THAT(actual_in_place, testing::ContainerEq(expected))
        << "From " << testing::PrintToString(lhs) << " | " << testing::PrintToString(rhs);
  });
}
TYPED_TEST_P(SetOpsTests, SetIntersection) {
  using elem_vec_t = std::vector<typename TypeParam::elem_t>;
  using container_t = typename TypeParam::container_t;
  using rhs_container_t = typename TypeParam::rhs_container_t;

  this->for_all_example_pairs([this](const elem_vec_t& lhs, const elem_vec_t& rhs) {
    container_t lhs_c = TypeParam::container_for_elems(lhs);
    rhs_container_t rhs_c = TypeParam::rhs_container_for_elems(rhs);
    elem_vec_t expected = this->set_intersection(lhs, rhs);

    elem_vec_t actual = TypeParam::elems_in_container(lhs_c & rhs_c);
    EXPECT_THAT(actual, testing::ContainerEq(expected))
        << "From " << testing::PrintToString(lhs) << " & " << testing::PrintToString(rhs);

    auto* in_place_result = &(lhs_c &= rhs_c);
    EXPECT_EQ(in_place_result, &lhs_c);
    elem_vec_t actual_in_place = TypeParam::elems_in_container(lhs_c);
    EXPECT_THAT(actual_in_place, testing::ContainerEq(expected))
        << "From " << testing::PrintToString(lhs) << " & " << testing::PrintToString(rhs);
  });
}

TYPED_TEST_P(SetOpsTests, SetDifference) {
  using elem_vec_t = std::vector<typename TypeParam::elem_t>;
  using container_t = typename TypeParam::container_t;
  using rhs_container_t = typename TypeParam::rhs_container_t;

  this->for_all_example_pairs([this](const elem_vec_t& lhs, const elem_vec_t& rhs) {
    container_t lhs_c = TypeParam::container_for_elems(lhs);
    rhs_container_t rhs_c = TypeParam::rhs_container_for_elems(rhs);
    elem_vec_t expected = this->set_difference(lhs, rhs);

    elem_vec_t actual = TypeParam::elems_in_container(lhs_c - rhs_c);
    EXPECT_THAT(actual, testing::ContainerEq(expected))
        << "From " << testing::PrintToString(lhs) << " - " << testing::PrintToString(rhs);

    auto* in_place_result = &(lhs_c -= rhs_c);
    EXPECT_EQ(in_place_result, &lhs_c);
    elem_vec_t actual_in_place = TypeParam::elems_in_container(lhs_c);
    EXPECT_THAT(actual_in_place, testing::ContainerEq(expected))
        << "From " << testing::PrintToString(lhs) << " - " << testing::PrintToString(rhs);
  });
}

REGISTER_TYPED_TEST_CASE_P(SetOpsTests, Empty, CreateAndAccess,  //
                           SetUnion, SetIntersection, SetDifference);

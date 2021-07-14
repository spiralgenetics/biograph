#include <set>

#include "modules/test/set_ops.h"

// Simple set bsaed on std::set that supports standard set operators.
// Implementations should subclass this, speciying type types of the
// sets that are on the left and right sides of the operators.
template <typename LhsType, typename RhsType>
struct my_set_tmpl : public std::set<int> {
 public:
  using std::set<int>::set;

  LhsType operator|(const RhsType& rhs) const {
    LhsType result;
    std::set_union(begin(), end(), rhs.begin(), rhs.end(), std::inserter(result, result.end()));
    return result;
  }
  LhsType operator&(const RhsType& rhs) const {
    LhsType result;
    std::set_intersection(begin(), end(), rhs.begin(), rhs.end(),
                          std::inserter(result, result.end()));
    return result;
  }
  LhsType operator-(const RhsType& rhs) const {
    LhsType result;
    std::set_difference(begin(), end(), rhs.begin(), rhs.end(),
                        std::inserter(result, result.end()));
    return result;
  }

  LhsType& operator|=(const RhsType& rhs) {
    (*this) = (*this) | rhs;
    return static_cast<LhsType&>(*this);
  }
  LhsType& operator&=(const RhsType& rhs) {
    (*this) = (*this) & rhs;
    return static_cast<LhsType&>(*this);
  }
  LhsType& operator-=(const RhsType& rhs) {
    (*this) = (*this) - rhs;
    return static_cast<LhsType&>(*this);
  }
};

// Set operations between set and itself.
struct my_set : my_set_tmpl<my_set, my_set> {
  using my_set_tmpl<my_set, my_set>::my_set_tmpl;
};

struct my_set_traits : set_ops_test_traits_base<my_set_traits, int, my_set> {
  static std::vector<int> example_elems() { return {1, 2, 3, 4}; };
};

INSTANTIATE_TYPED_TEST_CASE_P(MySet, SetOpsTests, my_set_traits);

// Set operations between set and a different type.
struct my_rhs_set : std::set<int> {
  using std::set<int>::set;
};

struct my_lhs_set : my_set_tmpl<my_lhs_set, my_rhs_set> {
  using my_set_tmpl<my_lhs_set, my_rhs_set>::my_set_tmpl;
};

struct my_lhs_rhs_set_traits
    : set_ops_test_traits_base<my_lhs_rhs_set_traits, int, my_lhs_set, my_rhs_set> {
  static std::vector<int> example_elems() { return {1, 2, 3, 4}; };
};

INSTANTIATE_TYPED_TEST_CASE_P(MyLhsRhsSet, SetOpsTests, my_lhs_rhs_set_traits);

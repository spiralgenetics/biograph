#include "modules/variants/sort.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class sort_test : public assemble_test,
                  public WithParamInterface<assembly::ordering_t> {
 public:
  void SetUp() { m_sort.emplace(GetParam(), test_output()); }
  void add(assembly a) { m_sort->add(make_unique<assembly>(a)); }
  void flush() {
    m_sort.reset();
    expect_sorted(GetParam());
  }

 protected:
  boost::optional<sorter> m_sort;
};

TEST_P(sort_test, sorts) {
  assembly a;
  a.right_offset = 100;
  a.seq = tseq("abc");

  a.left_offset = 5;
  add(a);
  a.left_offset = 3;
  add(a);
  a.left_offset = 50;
  add(a);

  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIs(3, tseq("abc"), 100),
                                        AssemblyIs(5, tseq("abc"), 100),
                                        AssemblyIs(50, tseq("abc"), 100)));
}

INSTANTIATE_TEST_CASE_P(left_offset_sort_test, sort_test,
                        ::testing::Values(assembly::left_offset_less_than));
INSTANTIATE_TEST_CASE_P(left_anchor_end_sort_test, sort_test,
                        ::testing::Values(assembly::left_anchor_end_less_than));

}  // namespace variants

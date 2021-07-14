#include "modules/variants/dedup.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;
class dedup_test : public assemble_test {
 public:
  void SetUp() { m_dedup.emplace(test_output()); }
  void add(assembly a) { m_dedup->add(make_unique<assembly>(a)); }
  void flush() {
    m_dedup.reset();
    expect_sorted(assembly::left_offset_less_than);
  }

 protected:
  boost::optional<deduper> m_dedup;
};

TEST_F(dedup_test, duplicate) {
  assembly a;
  a.left_offset = 1;
  a.right_offset = 5;
  a.seq = tseq("abc");
  add(a);
  add(a);

  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIs(1, tseq("abc"), 5)));
}

TEST_F(dedup_test, non_duplicate) {
  assembly a;
  a.left_offset = 1;
  a.right_offset = 5;
  a.seq = tseq("abc");
  add(a);

  // This one is duplicated.
  a.seq = tseq("def");
  add(a);
  add(a);

  a.right_offset = 6;
  add(a);
  a.left_offset = 2;
  add(a);

  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(AssemblyIs(1, tseq("abc"), 5),
                                   AssemblyIs(1, tseq("def"), 5),
                                   AssemblyIs(1, tseq("def"), 6),
                                   AssemblyIs(2, tseq("def"), 6)));
  ASSERT_TRUE(!m_assemblies.empty());
  // Make sure the one at the next position comes last.
  EXPECT_EQ(m_assemblies.back().left_offset, 2);
}

TEST_F(dedup_test, overlapping) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});

  {
    assembly a;
    a.left_offset = 0;
    a.left_anchor_len = tseq("ab").size();
    a.seq = tseq("ab") + dna_G + tseq("d");
    a.right_offset = tseq("abcd").size();
    a.right_anchor_len = tseq("d").size();
    add(a);
  }

  {
    assembly a;
    a.left_offset = tseq("a").size();
    a.left_anchor_len = tseq("b").size();
    a.seq = tseq("b") + dna_G + tseq("de");
    a.right_offset = tseq("abcde").size();
    a.right_anchor_len = tseq("de").size();
    add(a);
  }

  flush();
  ASSERT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(0, tseq("ab") + dna_G + tseq("de"),
                                     tseq("abcde").size())));
  const assembly& a0 = m_assemblies[0];
  EXPECT_EQ(a0.seq, tseq("ab") + dna_G + tseq("de"));
  EXPECT_EQ(a0.left_anchor_len, tseq("ab").size());
  EXPECT_EQ(a0.right_anchor_len, tseq("de").size());
}

TEST_F(dedup_test, encompassing) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});

  {
    assembly a;
    a.left_offset = 0;
    a.left_anchor_len = tseq("ab").size();
    a.seq = tseq("ab") + dna_G + tseq("de");
    a.right_offset = tseq("abcde").size();
    a.right_anchor_len = tseq("de").size();
    add(a);
  }

  {
    assembly a;
    a.left_offset = tseq("a").size();
    a.left_anchor_len = tseq("b").size();
    a.seq = tseq("b") + dna_G + tseq("d");
    a.right_offset = tseq("abcd").size();
    a.right_anchor_len = tseq("d").size();
    add(a);
  }

  flush();
  ASSERT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(0, tseq("ab") + dna_G + tseq("de"),
                                     tseq("abcde").size())));
  const assembly& a0 = m_assemblies[0];
  EXPECT_EQ(a0.seq, tseq("ab") + dna_G + tseq("de"));
  EXPECT_EQ(a0.left_anchor_len, tseq("ab").size());
  EXPECT_EQ(a0.right_anchor_len, tseq("de").size());
}

}  // namespace variants

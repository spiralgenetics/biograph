#include "modules/variants/trim_ref.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class trim_ref_test : public assemble_test {
 public:
  trim_ref_test() { m_trim_ref.emplace(m_options, test_output()); }
  void add(assembly a) { m_trim_ref->add(make_unique<assembly>(a)); }
  void flush() { m_trim_ref.reset(); }

 protected:
  boost::optional<ref_trimmer> m_trim_ref;
};

TEST_F(trim_ref_test, ignores_reference_passing_expand) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  assembly a;
  a.matches_reference = true;
  a.left_offset = 0;
  a.right_offset = tseq("abc").size();
  a.seq = tseq("abc");
  add(a);

  flush();

  EXPECT_THAT(m_assemblies, IsEmpty());
}

TEST_F(trim_ref_test, ignores_variant_equiv_to_reference) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  assembly a;
  a.matches_reference = false;
  a.left_offset = 0;
  a.right_offset = tseq("abc").size();
  a.seq = tseq("abc");
  add(a);

  flush();

  EXPECT_THAT(m_assemblies, IsEmpty());
}

TEST_F(trim_ref_test, expands) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});

  assembly a;
  a.left_offset = tseq("abc").size();
  a.seq = tseq("de") + dna_G + tseq("FGHI") + dna_G + tseq("jklmn");
  a.left_anchor_len = tseq("d").size();
  a.right_offset = tseq("abcdefghijklmn").size();
  a.right_anchor_len = tseq("mn").size();

  add(a);
  flush();

  ASSERT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(tseq("abcde").size(), dna_G + tseq("FGHI") + dna_G,
                                     tseq("abcdefghi").size())));
  const auto& a0 = m_assemblies[0];
  EXPECT_EQ(a0.left_anchor_len, 0);
  EXPECT_EQ(a0.right_anchor_len, 0);
}

TEST_F(trim_ref_test, duplicate_insert_expand) {
  // Duplicate middle should call as far left as possible.
  use_ref_parts({{0, tseq("abc") + dna_G + tseq("defghijklmnopqrstuvwxyz")}});

  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = 1;
  a.right_offset = (tseq("abc") + dna_G + tseq("defghi")).size();
  a.right_anchor_len = 1;
  a.seq = tseq("abc") + dna_G + dna_G + tseq("defghi");

  add(a);
  flush();

  ASSERT_THAT(m_assemblies, ElementsAre(AssemblyIs(tseq("abc").size(), dna_G, tseq("abc").size())));
  const assembly& a0 = m_assemblies[0];
  EXPECT_EQ(a0.right_anchor_len, 0);
  EXPECT_EQ(a0.left_anchor_len, 0);
}

TEST_F(trim_ref_test, duplicate_delete_expand) {
  // Duplicate middle should call as far left as possible.
  use_ref_parts({{0, tseq("abc") + dna_G + dna_G + tseq("defghijklmnopqrstuvwxyz")}});

  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = 1;
  a.right_offset = (tseq("abc") + dna_G + dna_G + tseq("defghi")).size();
  a.right_anchor_len = 1;
  a.seq = tseq("abc") + dna_G + tseq("defghi");

  add(a);
  flush();

  ASSERT_THAT(m_assemblies, ElementsAre(AssemblyIs(tseq("abc").size(), dna_sequence(),
                                                   (tseq("abc") + dna_G).size())));
  const assembly& a0 = m_assemblies[0];
  EXPECT_EQ(a0.left_anchor_len, 0);
  EXPECT_EQ(a0.right_anchor_len, 0);
}

TEST_F(trim_ref_test, left_anchor) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});

  assembly a;
  a.left_offset = tseq("abcd").size();
  a.right_offset = optional_aoffset::none;
  a.seq = tseq("e") + dna_G;

  add(a);
  flush();

  ASSERT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(tseq("abcde").size(), dna_G, optional_aoffset::none)));
}

TEST_F(trim_ref_test, right_anchor) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});

  assembly a;
  a.left_offset = optional_aoffset::none;
  a.right_offset = tseq("abcd").size();
  a.seq = dna_G + tseq("d");

  add(a);
  flush();

  ASSERT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(optional_aoffset::none, dna_G, tseq("abc").size())));
}

}  // namespace variants

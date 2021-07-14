#include "modules/variants/add_ref.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class add_ref_test : public assemble_test {
 public:
  void start_adder(aoffset_t pad_size, bool whole_ref = false, int max_len = 0) {
    m_adder.emplace(m_options, pad_size, whole_ref, max_len, test_output());
  }
  void add(aoffset_t left_offset, dna_sequence seq, aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    m_expected_non_ref.push_back(AssemblyIs(a->left_offset, a->seq, a->right_offset));
    m_adder->add(std::move(a));
  }

  void flush() {
    m_adder.reset();
    expect_sorted(assembly::left_offset_less_than);

    // add_ref should not add any non-reference assemblies.
    EXPECT_THAT(m_non_ref_assemblies, ElementsAreArray(m_expected_non_ref));
  }

 protected:
  boost::optional<add_ref> m_adder;
  std::vector<Matcher<const assembly&>> m_expected_non_ref;
};

TEST_F(add_ref_test, simple) {
  use_ref_parts({{0, dna_G + tseq("abcdefghijklmnopqrstuvwxyz")}});

  start_adder(0 /* pad size */);
  add(dna_G.size(), dna_T, (dna_G + tseq("a")).size());
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(dna_G.size(), (dna_G + tseq("a")).size())));
}

TEST_F(add_ref_test, pad_ref_ege_overlap) {
  use_ref_parts({{0, tseq("abcdefg")}});

  start_adder(tseq("a").size() + 1);
  add(tseq("a").size(), dna_T, tseq("abcdef").size());
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(0, tseq("a").size()),
                          RefAssemblyIs(tseq("a").size(), tseq("abcdef").size()),
                          RefAssemblyIs(tseq("abcdef").size(), tseq("abcdefg").size())));
}

TEST_F(add_ref_test, pad_ref_edge_underlap) {
  use_ref_parts({{0, tseq("abcdefg")}});

  start_adder(tseq("a").size());
  add(tseq("a").size(), dna_T, tseq("abcdef").size());
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(0, tseq("a").size()),
                          RefAssemblyIs(tseq("a").size(), tseq("abcdef").size()),
                          RefAssemblyIs(tseq("abcdef").size(), tseq("abcdefg").size())));
}

TEST_F(add_ref_test, pad_ref_edge_exact) {
  use_ref_parts({{0, dna_G + tseq("abcdefg") + dna_G}});

  start_adder(tseq("a").size());
  add(1 + tseq("a").size(), dna_T, 1 + tseq("abcdef").size());
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(1, 1 + tseq("a").size()),
                          RefAssemblyIs(1 + tseq("a").size(), 1 + tseq("abcdef").size()),
                          RefAssemblyIs(1 + tseq("abcdef").size(), 1 + tseq("abcdefg").size())));
}

TEST_F(add_ref_test, pad_asm_overlap) {
  use_ref_parts({{0, tseq("abc") + dna_G + tseq("def")}});

  start_adder(tseq("a").size() + 1);
  add(0, dna_T, tseq("ab").size());
  add(tseq("abc").size() + 1 + tseq("d").size(), dna_T,
      tseq("abc").size() + 1 + tseq("def").size());
  flush();

  EXPECT_THAT(
      m_ref_assemblies,
      ElementsAre(RefAssemblyIs(0, tseq("ab").size()),
                  RefAssemblyIs(tseq("ab").size(), tseq("abc").size() + 1 + tseq("d").size()),
                  RefAssemblyIs(tseq("abc").size() + 1 + tseq("d").size(),
                                tseq("abc").size() + 1 + tseq("def").size())));
}

TEST_F(add_ref_test, pad_asm_underlap) {
  use_ref_parts({{0, tseq("abc") + dna_G + tseq("def")}});

  start_adder(tseq("a").size());
  add(0, dna_T, tseq("ab").size());
  add(tseq("abc").size() + 1 + tseq("d").size(), dna_T,
      tseq("abc").size() + 1 + tseq("def").size());
  flush();

  EXPECT_THAT(
      m_ref_assemblies,
      ElementsAre(RefAssemblyIs(0, tseq("ab").size()),
                  RefAssemblyIs(tseq("ab").size(), tseq("abc").size()),
                  RefAssemblyIs(tseq("abc").size() + 1, tseq("abc").size() + 1 + tseq("d").size()),
                  RefAssemblyIs(tseq("abc").size() + 1 + tseq("d").size(),
                                tseq("abc").size() + 1 + tseq("def").size())));
}
TEST_F(add_ref_test, pad_asm_exact) {
  use_ref_parts({{0, tseq("abcdef")}});

  start_adder(tseq("a").size());
  add(0, dna_T, tseq("ab").size());
  add(tseq("abcd").size(), dna_T, tseq("abcdef").size());
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(0, tseq("ab").size()),
                          RefAssemblyIs(tseq("ab").size(), tseq("abcd").size()),
                          RefAssemblyIs(tseq("abcd").size(), tseq("abcdef").size())));
}

TEST_F(add_ref_test, overlapping_asm) {
  use_ref_parts({{0, tseq("abcdef")}});

  start_adder(tseq("a").size());
  add(0, dna_T, tseq("abcd").size());
  add(tseq("abc").size(), dna_T, tseq("abcdef").size());
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(0, tseq("abc").size()),
                          RefAssemblyIs(tseq("abc").size(), tseq("abcd").size()),
                          RefAssemblyIs(tseq("abcd").size(), tseq("abcdef").size())));
}

TEST_F(add_ref_test, whole_ref) {
  use_ref_parts({{0, tseq("abcdef")}});

  start_adder(tseq("a").size(), true /* whole ref */);
  flush();

  EXPECT_THAT(m_ref_assemblies, ElementsAre(RefAssemblyIs(0, tseq("abcdef").size())));
}

TEST_F(add_ref_test, whole_ref_multiple_scaffolds) {
  use_ref_parts({{10, tseq("abc")}, {1000, tseq("def")}});

  start_adder(0, true /* whole ref */);
  flush();

  EXPECT_THAT(m_ref_assemblies, ElementsAre(RefAssemblyIs(10, 10 + tseq("abc").size()),
                                            RefAssemblyIs(1000, 1000 + tseq("def").size())));
}

TEST_F(add_ref_test, whole_ref_with_variant) {
  use_ref_parts({{0, tseq("abcdef")}});

  start_adder(1, true /* whole ref */);
  add(tseq("ab").size(), dna_T, tseq("abc").size());
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(0, tseq("ab").size()),
                          RefAssemblyIs(tseq("ab").size(), tseq("abc").size()),
                          RefAssemblyIs(tseq("abc").size(), tseq("abcdef").size())));
}

TEST_F(add_ref_test, whole_ref_max_len) {
  use_ref_parts({{0, tseq("abcdef")}});

  start_adder(0, true /* whole ref */, tseq("abcd").size() /* max len */);
  flush();

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(0, tseq("abcd").size()),
                          RefAssemblyIs(tseq("abcd").size(), tseq("abcdef").size())));
}

}  // namespace variants

#include "modules/variants/assemble.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

using namespace testing;
using namespace variants;
using namespace dna_testutil;

class merge_test : public Test {
 public:
  aligned_var av(aoffset_t left_offset, aoffset_t right_offset, dna_sequence seq) {
    aligned_var result;
    result.left_offset = left_offset;
    result.right_offset = right_offset;
    result.seq = seq;
    return result;
  }
  assembly_ptr make_assembly(aoffset_t left_offset, aoffset_t right_offset, dna_sequence seq,
                             std::vector<aligned_var> aligned_vars = {}) {
    assembly_ptr a = make_unique<assembly>();
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    a->aligned_variants = aligned_vars;
    if (aligned_vars.empty()) {
      a->matches_reference = true;
    }
    return a;
  }

  void do_merge() { m_result = merge_assemblies(*m_in1, *m_in2); }

 protected:
  assembly_ptr m_in1;
  assembly_ptr m_in2;

  assembly_ptr m_result;
};

TEST_F(merge_test, disjoint) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(41, 71, tseq("ABC"));
  do_merge();
  EXPECT_FALSE(m_result) << *m_result;
}

TEST_F(merge_test, touching) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(40, 70, tseq("DEF"));
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, RefAssemblyIs(10, 70));
  EXPECT_EQ(m_result->seq, tseq("ABCDEF"));
}

TEST_F(merge_test, bad_overlap) {
  // This normally wouldn't happen since they'd both be matching the same reference
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(30, 60, tseq("DEF"));
  do_merge();
  EXPECT_FALSE(m_result) << *m_result;
}

TEST_F(merge_test, good_overlap) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(30, 60, tseq("CDE"));
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, RefAssemblyIs(10, 60));
  EXPECT_EQ(m_result->seq, tseq("ABCDE"));
}

TEST_F(merge_test, whole_overlap) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(10, 40, tseq("ABC"));
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, RefAssemblyIs(10, 40));
  EXPECT_EQ(m_result->seq, tseq("ABC"));
}

TEST_F(merge_test, whole_overlap_bad) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(10, 40, tseq("AbbC"), {av(20, 30, tseq("bb"))});
  do_merge();
  EXPECT_FALSE(m_result);
}

TEST_F(merge_test, subsumed_bad) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(20, 30, tseq("D"));
  do_merge();
  EXPECT_FALSE(m_result) << *m_result;
}

TEST_F(merge_test, subsumed_good) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(20, 30, tseq("B"));
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, RefAssemblyIs(10, 40));
  EXPECT_EQ(m_result->seq, tseq("ABC"));
}

TEST_F(merge_test, disjoint_vars_good) {
  m_in1 = make_assembly(10, 40, tseq("AbbC"), {av(20, 30, tseq("bb"))});
  m_in2 = make_assembly(40, 70, tseq("DeeF"), {av(50, 60, tseq("ee"))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("AbbCDeeF"), 70));
  EXPECT_THAT(m_result->aligned_variants,
              ElementsAre(av(20, 30, tseq("bb")), av(50, 60, tseq("ee"))));
}

TEST_F(merge_test, overlapping_vars_good) {
  m_in1 = make_assembly(10, 50, tseq("ABccD"), {av(30, 40, tseq("cc"))});
  m_in2 = make_assembly(20, 60, tseq("BccDE"), {av(30, 40, tseq("cc"))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("ABccDE"), 60));
  EXPECT_THAT(m_result->aligned_variants, ElementsAre(av(30, 40, tseq("cc"))));
}

TEST_F(merge_test, overlapping_vars_bad) {
  m_in1 = make_assembly(10, 50, tseq("ABccD"), {av(30, 40, tseq("cc"))});
  m_in2 = make_assembly(20, 60, tseq("BxxDE"), {av(30, 40, tseq("xx"))});
  do_merge();
  EXPECT_FALSE(m_result) << *m_result;
}

TEST_F(merge_test, subsumed_vars_good) {
  m_in1 = make_assembly(10, 80, tseq("AbbCddEffG"),
                        {av(20, 30, tseq("bb")), av(40, 50, tseq("dd")), av(60, 70, tseq("ff"))});
  m_in2 = make_assembly(30, 60, tseq("CddE"), {av(40, 50, tseq("dd"))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("AbbCddEffG"), 80));
  EXPECT_THAT(m_result->aligned_variants,
              ElementsAre(av(20, 30, tseq("bb")), av(40, 50, tseq("dd")), av(60, 70, tseq("ff"))));
}

TEST_F(merge_test, left_ref) {
  m_in1 = make_assembly(10, 40, tseq("ABC"));
  m_in2 = make_assembly(30, 60, tseq("CddE"), {av(40, 50, tseq("dd"))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("ABCddE"), 60));
  EXPECT_THAT(m_result->aligned_variants, ElementsAre(av(40, 50, tseq("dd"))));
}

TEST_F(merge_test, right_ref) {
  m_in1 = make_assembly(10, 40, tseq("AbbC"), {av(20, 30, tseq("bb"))});
  m_in2 = make_assembly(30, 60, tseq("CDE"));
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("AbbCDE"), 60));
  EXPECT_THAT(m_result->aligned_variants, ElementsAre(av(20, 30, tseq("bb"))));
}

TEST_F(merge_test, delete_overlap) {
  m_in1 = make_assembly(10, 70, tseq("ABF"), {av(30, 60, tseq(""))});
  m_in2 = make_assembly(20, 80, tseq("BFG"), {av(30, 60, tseq(""))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("ABFG"), 80));
  EXPECT_THAT(m_result->aligned_variants, ElementsAre(av(30, 60, tseq(""))));
}

TEST_F(merge_test, delete_adjacent) {
  m_in1 = make_assembly(10, 100, tseq("A"), {av(20, 100, tseq(""))});
  m_in2 = make_assembly(100, 200, tseq("H"), {av(100, 190, tseq(""))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("AH"), 200));
  EXPECT_THAT(m_result->aligned_variants,
              ElementsAre(av(20, 100, tseq("")), av(100, 190, tseq(""))));
}

TEST_F(merge_test, insert_right_adjacent) {
  m_in1 = make_assembly(10, 100, tseq("A"), {av(20, 100, tseq(""))});
  m_in2 = make_assembly(100, 110, tseq("efgH"), {av(100, 100, tseq("efg"))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(10, tseq("AefgH"), 110));
  EXPECT_THAT(m_result->aligned_variants,
              ElementsAre(av(20, 100, tseq("")), av(100, 100, tseq("efg"))));
}

TEST_F(merge_test, insert_both_adjacent) {
  m_in1 = make_assembly(90, 100, tseq("Abcd"), {av(100, 100, tseq("bcd"))});
  m_in2 = make_assembly(100, 110, tseq("efgH"), {av(100, 100, tseq("efg"))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(90, tseq("AbcdefgH"), 110));
  EXPECT_THAT(m_result->aligned_variants,
              ElementsAre(av(100, 100, tseq("bcd")), av(100, 100, tseq("efg"))));
}

TEST_F(merge_test, insert_left_adjacent) {
  m_in1 = make_assembly(90, 100, tseq("Abcd"), {av(100, 100, tseq("bcd"))});
  m_in2 = make_assembly(100, 200, tseq("H"), {av(100, 190, tseq(""))});
  do_merge();
  ASSERT_TRUE(m_result);
  EXPECT_THAT(*m_result, AssemblyIs(90, tseq("AbcdH"), 200));
  EXPECT_THAT(m_result->aligned_variants,
              ElementsAre(av(100, 100, tseq("bcd")), av(100, 190, tseq(""))));
}

// Merges seen in the wild
TEST_F(merge_test, wild1) {
  aligned_var v;

  m_in1 = make_unique<assembly>();
  m_in1->left_offset = 248033565;
  m_in1->right_offset = 248033882;
  m_in1->seq = dna_sequence(
      "CATGCCCAGAAAAATTAGGATTGGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCATAGAGCCTGCACCACCATGCCCAGAAGAACGA"
      "GGATTAGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCATAGAGCCTGCACCACCATGCCCAGAAGAACGAGGATTAGTATATTCATG"
      "AGCAGTATCAGAGTCACGGCCACCCATAGAGCCTGCACCACCACACCCAGAAGAACGAGGATTAGTATATTCATGAGCAGTATCAGAGTCAC"
      "GGCCACCCGTAGAGCCTGCACCACCACGCCCAGAAGAACGAGGATTAGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCATAGAGCCT"
      "GCACCACCACGCCCAGAAGAATGAGGATTAGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCATAGAGCCTGCACCACCACACCCAGA"
      "AGAACGAGGATTAGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCATAGAGCCTGCACCACCACGCCCAGAAGAACGA");

  v.left_offset = 248033717;
  v.right_offset = 248033718;
  v.seq = dna_sequence("T");
  m_in1->aligned_variants.push_back(v);

  v.left_offset = 248033730;
  v.right_offset = 248033730;
  v.seq = dna_sequence(
      "GAGGATTAGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCATAGAGCCTGCACCACCACACCCAGAAGAACGAGGATTAGTATATTCA"
      "TGAGCAGTATCAGAGTCACGGCCACCCGTAGAGCCTGCACCACCACGCCCAGAAGAACGAGGATTAGTATATTCATGAGCAGTATCAGAGTC"
      "ACGGCCACCCATAGAGCCTGCACCACCACGCCCAGAAGAAT");
  m_in1->aligned_variants.push_back(v);

  m_in2 = make_unique<assembly>();
  m_in2->left_offset = 248033766;
  m_in2->right_offset = 248034153;
  m_in2->seq = dna_sequence(
      "GGCCACCCATAGAGCCTGCACCACCACACCCAGAAGAACGAGGATTAGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCATAGAGCCT"
      "GCACCACCACGCCCAGAAGAACGAGGATTAGTATATTCATGAGCAGTATCAGAGTCACGGCCACCCGTAGAGCCTGCACCACCACGCCCAGA"
      "AGAACGAGGATTAGTATATTCATGAGCAGTATCAGAGTCATGGCCACCCATAG");

  v.left_offset = 248033794;
  v.right_offset = 248033869;
  v.seq = dna_sequence("");
  m_in2->aligned_variants.push_back(v);

  v.left_offset = 248033925;
  v.right_offset = 248034000;
  v.seq = dna_sequence("");
  m_in2->aligned_variants.push_back(v);

  do_merge();
  EXPECT_FALSE(m_result) << *m_result;
}

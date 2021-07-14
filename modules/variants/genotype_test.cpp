#include "modules/variants/genotype.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

// TODO(nils): Figure out what to do and add a test to do teh right
// thing in the case where we have a deletion that spans a section of
// reference that includes an ambiguous portion that has coverage
// elsewhere in the reference that eclipses that of the deletion.
namespace variants {

using namespace testing;
using namespace dna_testutil;
using namespace coverage_testutil;

class genotype_test : public assemble_test {
 public:
  aligned_var av(aoffset_t left_offset, aoffset_t right_offset, dna_sequence seq) {
    aligned_var result;
    result.left_offset = left_offset;
    result.right_offset = right_offset;
    result.seq = seq;
    return result;
  }
  assembly_ptr make_as(aoffset_t left_offset, aoffset_t right_offset,
                       std::vector<aligned_var> aligned_vars, coverage_constructor cov) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = ++m_assembly_id;
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->aligned_variants = aligned_vars;
    a->coverage = cov;
    a->matches_reference = aligned_vars.empty();

    // Reconstruct seq from reference and vars
    aoffset_t ref_offset = a->left_offset;
    auto vit = aligned_vars.begin();
    while (vit != aligned_vars.end()) {
      aoffset_t ref_adv = vit->left_offset - ref_offset;
      a->seq += m_scaffold.subscaffold(ref_offset, ref_adv).get_simple();
      a->seq += vit->seq;
      ref_offset = vit->right_offset;
      ++vit;
    }
    aoffset_t ref_adv = a->right_offset - ref_offset;
    a->seq += m_scaffold.subscaffold(ref_offset, ref_adv).get_simple();

    return a;
  }
  void add_as(aoffset_t left_offset, aoffset_t right_offset, std::vector<aligned_var> aligned_vars,
              coverage_constructor cov) {
    add(make_as(left_offset, right_offset, aligned_vars, cov));
  }

  void add(assembly_ptr a) {
    CHECK(m_genotyper);
    check_assembly(*a, "genotype_test");
    m_genotyper->add(std::move(a));
  }

  void start_genotyper() { m_genotyper = make_unique<genotyper>(m_options, test_output()); }

  void flush() {
    CHECK(m_genotyper);
    m_genotyper.reset();
  }

 protected:
  std::unique_ptr<genotyper> m_genotyper;
  size_t m_assembly_id = 0;
};

TEST_F(genotype_test, simple_ref) {
  use_ref_parts({{0, tseq("abcdefg")}});

  start_genotyper();
  add_as(0, tseq("abcdefg").size(), {}, 10 + over("abcdefg", 10) + 10);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(RefAssemblyIs(0, tseq("abcdefg").size()),
                                              CoverageIs(10 + over("abcdefg", 10) + 10))));
}

TEST_F(genotype_test, partial_ref) {
  use_ref_parts({{0, tseq("abcdefgh")}});

  start_genotyper();
  add_as(tseq("a").size(), tseq("abcdefg").size(), {},
         0 + over("b", 0) + 10 + over("cd", 10) + 10 + over("efg", 0) + 0);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(RefAssemblyIs(tseq("ab").size(), tseq("abcd").size()),
                                              CoverageIs(10 + over("cd", 10) + 10))));
}

TEST_F(genotype_test, single_insert) {
  use_ref_parts({{0, tseq("abcdefgh")}});

  start_genotyper();
  add_as(30, 50, {av(40, 40, dna_A)}, 0 + over("d", 0) +    //
                                          10 +              //
                                          /* inserted A */  //
                                          10 + over("e", 0) + 0);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(AssemblyIs(40, dna_A, 40), CoverageIs(rpt(2, 10)))));
}

TEST_F(genotype_test, two_inserts) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efgh")}});

  start_genotyper();
  add_as(30, 51, {av(40, 40, dna_A)}, 0 + over("d", 0) +      //
                                          10                  //
                                          /* inserted A */ +  //
                                          10                  //
                                          /* ref T */ +
                                          0 + over("e", 0) + 0);
  add_as(30, 51, {av(41, 41, dna_A)}, 0 + over("d", 0) +  //
                                          0               //
                                          /* ref T */ +
                                          10  //
                                          /* inserted A */ +
                                          10 + over("e", 0) + 0);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(AssemblyIs(40, dna_A, 40), CoverageIs(rpt(2, 10))),
                                        AllOf(AssemblyIs(41, dna_A, 41), CoverageIs(rpt(2, 10)))));
}

TEST_F(genotype_test, single_delete) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efgh")}});

  start_genotyper();
  add_as(30, 51, {av(40, 41, dna_sequence())}, 0 + over("d", 0) +  //
                                                   10 +            // Deletion
                                                   over("e", 0) + 0);
  flush();

  EXPECT_THAT(m_assemblies,
              ElementsAre(AllOf(AssemblyIs(40, dna_sequence(), 41), CoverageIs(rpt(1, 10)))));
}

TEST_F(genotype_test, abort_in_middle) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efgh")}});

  start_genotyper();
  add_as(30, 51, {av(40, 40, dna_A)}, 0 + over("d", 0) +      //
                                          10                  //
                                          /* inserted A */ +  //
                                          10                  //
                                          /* ref T */ +
                                          0 + over("e", 0) + 0);
  add_as(30, 51, {av(41, 41, dna_A)}, 0 + over("d", 0) +  //
                                          0               //
                                          /* ref T */ +
                                          10  //
                                          /* inserted A */ +
                                          10 + over("e", 0) + 0);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(AssemblyIs(40, dna_A, 40), CoverageIs(rpt(2, 10))),
                                        AllOf(AssemblyIs(41, dna_A, 41), CoverageIs(rpt(2, 10)))));
}

struct genotype_test_params {
  int pad_a_before = 0;
  int pad_a_after = 0;

  int pad_b_before = 0;
  int pad_b_after = 0;

  int b_offset = 0;
  int b_delete_size = 0;

  bool add_a_first = false;
  bool a_is_better = false;

  friend std::ostream& operator<<(std::ostream& os, const genotype_test_params& p) {
    return os << "A pad: " << p.pad_a_before << "," << p.pad_a_after
              << ", B pad: " << p.pad_b_before << "," << p.pad_b_after
              << ", B offset: " << p.b_offset << " B delete size: " << p.b_delete_size
              << ", adding " << (p.add_a_first ? "A" : "B") << " first, "
              << (p.a_is_better ? "A" : "B") << " is better.";
  }
};

class genotype_p_test : public genotype_test, public WithParamInterface<genotype_test_params> {
  void SetUp() override { m_p = GetParam(); }

 protected:
  genotype_test_params m_p;
};

// Tests conflicts between two assemblies, A and B, with a max ploid value of 1.
//
// A is a single insert at position 40
// B is a multiple base deletion of length b_delete_size.
TEST_P(genotype_p_test, ins_delete_conflict) {
  std::cout << "Running with params: " << m_p << "\n";
  std::cout.flush();
  int a_depth = m_p.a_is_better ? 10 : 5;
  int b_depth = m_p.a_is_better ? 5 : 10;
  m_options.max_ploids = 1;
  use_ref_parts({{0, tseq("abcdefgh")}});

  aoffset_t a_var_start = 40;
  aoffset_t a_var_end = 40;
  assembly_ptr a = make_as(40 - m_p.pad_a_before, 40 + m_p.pad_a_after, {av(40, 40, dna_A)},
                           rpt(m_p.pad_a_before, 0) + rpt(2, a_depth) + rpt(m_p.pad_a_after, 0));

  aoffset_t b_var_start = 40 + m_p.b_offset;
  aoffset_t b_var_end = b_var_start + m_p.b_delete_size;
  assembly_ptr b = make_as(b_var_start - m_p.pad_b_before, b_var_end + m_p.pad_b_after,
                           {av(b_var_start, b_var_end, dna_sequence())},
                           rpt(m_p.pad_b_before, 0) + rpt(1, b_depth) + rpt(m_p.pad_b_after, 0));

  // Skip adding them out of order (but if they start the same place, test adding in both orders)
  if (m_p.add_a_first && a->left_offset > b->left_offset) {
    return;
  }
  if (!m_p.add_a_first && a->left_offset < b->left_offset) {
    return;
  }

  start_genotyper();
  if (m_p.add_a_first) {
    add(std::move(a));
    add(std::move(b));
  } else {
    add(std::move(b));
    add(std::move(a));
  }
  flush();

  Matcher<assembly> a_matcher =
      AllOf(AssemblyIs(a_var_start, dna_A, a_var_end), CoverageIs(rpt(2, a_depth)));
  Matcher<assembly> b_matcher =
      AllOf(AssemblyIs(b_var_start, dna_sequence(), b_var_end), CoverageIs(rpt(1, b_depth)));

  if (a_var_end < b_var_start) {
    // Disjoint, A first
    EXPECT_THAT(m_assemblies, ElementsAre(a_matcher, b_matcher));
  } else if (a_var_start > b_var_end) {
    // Disjoint, a second
    EXPECT_THAT(m_assemblies, ElementsAre(b_matcher, a_matcher));
  } else if (m_p.a_is_better) {
    // Conflicts, A is better
    EXPECT_THAT(m_assemblies, ElementsAre(a_matcher));
  } else {
    // Conflicts, B is better
    EXPECT_THAT(m_assemblies, ElementsAre(b_matcher));
  }
}

// Full cartesian product of all test parameters, to try to get all
// the edge cases.
std::vector<genotype_test_params> all_params() {
  genotype_test_params p;
  std::vector<genotype_test_params> result;

#define FORALL(VARNAME, ...)         \
  for (auto VARNAME : __VA_ARGS__) { \
    p.VARNAME = VARNAME;
#define ENDFORALL }
  FORALL(pad_a_before, {0, 1, 2});
  FORALL(pad_b_before, {0, 1, 2});
  FORALL(pad_a_after, {0, 1, 2});
  FORALL(pad_b_after, {0, 1, 2});
  FORALL(b_delete_size, {1, 2, 3});
  for (p.b_offset = -p.b_delete_size - 2; p.b_offset <= p.b_delete_size + 2; ++p.b_offset) {
    FORALL(add_a_first, {false, true});
    FORALL(a_is_better, {false, true});
    result.push_back(p);
    ENDFORALL;
    ENDFORALL;
  }
  ENDFORALL;
  ENDFORALL;
  ENDFORALL;
  ENDFORALL;
  ENDFORALL;
#undef FORALL
#undef ENDFORALL
  return result;
}

INSTANTIATE_TEST_CASE_P(genotype_p_tests, genotype_p_test, ::testing::ValuesIn(all_params()));

}  // namespace variants

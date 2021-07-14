#include "modules/variants/split_variants.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/variants/assemble_testutil.h"

using namespace testing;
using namespace variants;
using namespace dna_testutil;

class split_variants_test : public assemble_test {
 public:
  aligned_var av(aoffset_t left_offset, aoffset_t right_offset, dna_sequence seq) {
    aligned_var result;
    result.left_offset = left_offset;
    result.right_offset = right_offset;
    result.seq = seq;
    return result;
  }
  assembly_ptr make_as(aoffset_t left_offset, aoffset_t right_offset,
                       std::vector<aligned_var> aligned_vars) {
    CHECK(m_split_variants);
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = ++m_assembly_id;
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->aligned_variants = aligned_vars;

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
  void add_as(aoffset_t left_offset, aoffset_t right_offset,
              std::vector<aligned_var> aligned_vars) {
    add(make_as(left_offset, right_offset, aligned_vars));
  }

  void add(assembly_ptr a) { m_split_variants->add(std::move(a)); }

  void flush() {
    m_split_variants.reset();

    // Make sure intervals either overlap completely or not at all.
    std::map<aoffset_t /* left offset */, int /* ref count */> ref_counts;
    for (const auto& a : m_assemblies) {
      ref_counts[a.left_offset] += a.matches_reference ? 1 : 0;
      EXPECT_THAT(a.aligned_variants, IsEmpty()) << a;

      for (const auto& b : m_assemblies) {
        if (b.left_offset < a.left_offset) {
          EXPECT_LE(b.right_offset, a.left_offset) << "a: " << a << " b: " << b;
        } else if (b.right_offset > a.right_offset) {
          EXPECT_GE(b.left_offset, a.right_offset) << "a: " << a << " b: " << b;
        } else {
          EXPECT_EQ(a.left_offset, b.left_offset) << "a: " << a << " b: " << b;
          EXPECT_EQ(a.right_offset, b.right_offset) << "a: " << a << " b: " << b;
        }
      }
    }
    EXPECT_THAT(ref_counts, Each(Pair(_, 1)));
  }

  void start_splitter() {
    // TODO(nils): Remove this option and rework the test not to need it.
    m_options.trace_reference_assemblies = true;

    m_split_variants = make_unique<split_variants>(m_options, test_output());
  }

 protected:
  std::unique_ptr<split_variants> m_split_variants;
  size_t m_assembly_id = 0;
};

template <typename T>
Matcher<const assembly&> OtherDepthIs(T&& t) {
  return Field(&assembly::other_depth, std::forward<T>(t));
}

TEST_F(split_variants_test, hetrozygous_snp) {
  use_ref_parts({{0, tseq("abcd") + dna_A + tseq("efghijklmnopqrstuvwxyz")}});
  // Reference coverage
  use_reads({
      // ref:
      tseq("abcd") + dna_A + tseq("efgh"),
      // Var:
      tseq("cd") + dna_T + tseq("efgh"), tseq("d") + dna_T + tseq("efghi"),
  });
  start_splitter();
  add_as(10, 100, {av(40, 41, dna_T)});
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(  //
                  AllOf(AssemblyIs(40, dna_T, 41), OtherDepthIs(2)),
                  AllOf(RefAssemblyIs(40, 41), OtherDepthIs(1))));
}

TEST_F(split_variants_test, homozygous_snp) {
  use_ref_parts({{0, tseq("abcd") + dna_A + tseq("efghijklmnopqrstuvwxyz")}});
  // Reference coverage
  use_reads({                                    // Var:
             tseq("cd") + dna_T + tseq("efgh"),  //
             dna_T + tseq("efghi"),
             tseq("cd") + dna_T + tseq("efgh"),  //
             tseq("d").subseq(9, 1) + dna_T + tseq("efghi"),
             // Var, but shouldn't be included since it doesn't span both sides:
             tseq("bcd").subseq(1, 19), dna_T + tseq("efghij"), tseq("efghijk")});
  start_splitter();
  add_as(10, 200, {av(40, 41, dna_T)});
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(  //
                  AllOf(AssemblyIs(40, dna_T, 41), OtherDepthIs(3)),
                  AllOf(RefAssemblyIs(40, 41), OtherDepthIs(0))));
}

TEST_F(split_variants_test, compound_hetrozygous_pad) {
  use_ref_parts({{0, tseq("abcd") + dna_A + dna_A + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({// Var 1
             tseq("cd") + dna_T + dna_A + tseq("e")});

  start_splitter();
  add_as(10, 100, {av(40, 41, dna_T)});
  add_as(20, 190, {av(41, 42, dna_T)});
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(  //
                  AllOf(AssemblyIs(40, dna_T + dna_A, 42), OtherDepthIs(1)),
                  AllOf(AssemblyIs(40, dna_A + dna_T, 42), OtherDepthIs(0)),
                  AllOf(RefAssemblyIs(40, 42), OtherDepthIs(0))));
}

TEST_F(split_variants_test, separate_inserts) {
  use_ref_parts({{0, tseq("abcd") + dna_A + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({// Var 1
             tseq("cd") + dna_T + dna_A + /* dna_T + */ tseq("e")});

  start_splitter();
  add_as(10, 100, {av(40, 40, dna_T)});
  add_as(20, 190, {av(41, 41, dna_T)});
  flush();

  EXPECT_THAT(m_assemblies,  //
              UnorderedElementsAre(
                  // First insert:
                  AllOf(AssemblyIs(39, dna_C + dna_T + dna_A, 41), OtherDepthIs(1)),
                  // Second insert:
                  AllOf(AssemblyIs(39, dna_C + dna_A + dna_T, 41), OtherDepthIs(0)),
                  AllOf(RefAssemblyIs(39, 41), OtherDepthIs(0))));
}

TEST_F(split_variants_test, multi_snps) {
  use_ref_parts(
      {{0, tseq("abcd") + dna_A + dna_A + dna_A + dna_A + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({// Var 1
             tseq("cd") + dna_T + dna_A + dna_A + dna_T + tseq("e")});

  start_splitter();
  add_as(10, 100, {// Separate:
                   av(40, 41, dna_T),
                   // Combined:
                   av(43, 44, dna_T)});
  // Combined:
  add_as(20, 190, {av(42, 43, dna_T)});
  flush();

  EXPECT_THAT(m_assemblies,  //
              UnorderedElementsAre(
                  // Separate snp:
                  AllOf(AssemblyIs(40, dna_T, 41), OtherDepthIs(1)),
                  // TODO(nils): Should we report a variant for this assembly even when it's the
                  // same as the reference?
                  AllOf(AssemblyIs(40, dna_A, 41), OtherDepthIs(0)),
                  AllOf(RefAssemblyIs(40, 41), OtherDepthIs(0)),
                  // Combined snps:
                  AllOf(AssemblyIs(42, dna_T + dna_A, 44), OtherDepthIs(0)),
                  AllOf(AssemblyIs(42, dna_A + dna_T, 44), OtherDepthIs(1)),
                  AllOf(RefAssemblyIs(42, 44), OtherDepthIs(0))));
}

TEST_F(split_variants_test, inserts) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({// Var 1
             tseq("cd") + dna_T + tseq("e")});

  start_splitter();
  add_as(10, 100, {av(40, 40, dna_T)});
  add_as(20, 190, {av(40, 40, dna_A)});
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(  //
                  AllOf(AssemblyIs(39, dna_C + dna_T, 40), OtherDepthIs(1)),
                  AllOf(AssemblyIs(39, dna_C + dna_A, 40), OtherDepthIs(0)),
                  AllOf(RefAssemblyIs(39, 40), OtherDepthIs(0))));
}

TEST_F(split_variants_test, close_inserts) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({// Var 1
             tseq("cd") + dna_A + dna_T + dna_A + tseq("e")});

  start_splitter();
  add_as(10, 100, {av(40, 40, dna_A), av(41, 41, dna_A)});
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(  //
                  AllOf(AssemblyIs(39, dna_C + dna_A + dna_T + dna_A, 41), OtherDepthIs(1)),
                  AllOf(RefAssemblyIs(39, 41), OtherDepthIs(0))));
}

TEST_F(split_variants_test, deletes) {
  use_ref_parts({{0, tseq("abcd") + dna_C + dna_T + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({// Var 1
             tseq("cd") + /* delete dna_C + */ dna_T + tseq("e")});

  start_splitter();
  add_as(10, 100, {av(40, 41, dna_sequence())});
  add_as(20, 190, {av(41, 42, dna_sequence())});
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(  //
                  AllOf(AssemblyIs(40, dna_T, 42), OtherDepthIs(1)),
                  AllOf(AssemblyIs(40, dna_C, 42), OtherDepthIs(0)),
                  AllOf(RefAssemblyIs(40, 42), OtherDepthIs(0))));
}

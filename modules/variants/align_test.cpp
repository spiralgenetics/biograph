#include "modules/variants/align.h"
#include "modules/variants/anchor_drop.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

MATCHER_P3(VarIs, start_offset, seq, limit_offset, "") {
  return arg.left_offset == aoffset_t(start_offset) &&
         arg.right_offset == aoffset_t(limit_offset) && arg.seq == seq;
}

template <typename... T>
Matcher<const assembly&> AlignedVarsAre(T&&... args) {
  return Field(&assembly::aligned_variants,
               ElementsAre(std::forward<T>(args)...));
}

template <typename T>
Matcher<const assembly&> AlignedVarsAreArray(T&& arg) {
  return Field(&assembly::aligned_variants,
               ElementsAreArray(std::forward<T>(arg)));
}

class align_test : public assemble_test {
 public:
  void align(assembly a) { align_multiple({a}); }
  void align_multiple(const std::vector<assembly> asms) {
    {
      auto splitter = make_unique<align_splitter>(test_output());
      auto save_aligned = make_unique<assemble_lambda_copy>(
          [this](const assembly& a) {
            m_aligned.push_back(a);
            std::cout << "Aligned: " << a << " with "
                      << a.aligned_variants.size() << " vars:\n";

            for (const auto& v : a.aligned_variants) {
              std::cout << "  " << v << "\n";
            }
            std::cout.flush();
          },
          std::move(splitter), "save_aligned");
      auto al = make_unique<aligner>(m_options, std::move(save_aligned));
      auto ad = make_unique<anchor_dropper>(m_options, std::move(al));
      for (const auto& a : asms) {
        ad->add(make_unique<assembly>(a));
      }
    }
    expect_sorted(assembly::left_offset_less_than);
  }

 protected:
  std::vector<assembly> m_aligned;
};

TEST_F(align_test, simple) {
  dna_sequence orig = dna_T;
  dna_sequence variant = dna_A + dna_A;

  m_scaffold.add(0, tseq("abc") + orig + tseq("def"));
  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = tseq("ab").size();
  a.right_offset = (tseq("abc") + orig + tseq("def")).size();
  a.right_anchor_len = tseq("ef").size();
  a.seq = tseq("abc") + variant + tseq("def");
  align(a);

  EXPECT_THAT(m_aligned,
              ElementsAre(AlignedVarsAre(VarIs(tseq("abc").size(), variant,
                                               (tseq("abc") + orig).size()))));

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(
          RefAssemblyIs(0, tseq("abc").size()),
          AssemblyIs(tseq("abc").size(), variant, (tseq("abc") + orig).size()),
          RefAssemblyIs((tseq("abc") + orig).size(),
                        (tseq("abc") + orig + tseq("def")).size())));
}

class anchor_drop_test : public align_test,
                         public WithParamInterface<unsigned /* kmer size */> {
 public:
  anchor_drop_test() {
    m_old_max_kmer_size =
        anchor_dropper::set_max_kmer_size_for_testing(GetParam());
  }
  ~anchor_drop_test() {
    CHECK_EQ(GetParam(), anchor_dropper::set_max_kmer_size_for_testing(
                             m_old_max_kmer_size));
  }

 private:
  unsigned m_old_max_kmer_size;
};

TEST_P(anchor_drop_test, drop_right) {
  dna_sequence orig = dna_T;
  dna_sequence variant = dna_A + dna_A;

  m_options.min_anchor_drop_overlap = tseq("d").size();
  m_scaffold.add(0, tseq("abc") + orig + tseq("defdefdef"));
  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = tseq("ab").size();
  a.right_offset = (tseq("abc") + orig + tseq("defdefdef")).size();
  a.right_anchor_len = 0;
  a.seq = tseq("abc") + variant + tseq("d");
  align(a);

  ASSERT_THAT(m_aligned,
              ElementsAre(AlignedVarsAre(VarIs(tseq("abc").size(), variant,
                                               (tseq("abc") + orig).size()))));
  EXPECT_EQ(m_aligned[0].left_offset, 0);
  EXPECT_EQ(m_aligned[0].right_offset, (tseq("abc") + orig + tseq("d")).size());

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(
          RefAssemblyIs(0, tseq("abc").size()),
          AssemblyIs(tseq("abc").size(), variant, (tseq("abc") + orig).size()),
          RefAssemblyIs((tseq("abc") + orig).size(),
                        (tseq("abc") + orig + tseq("d")).size())));
}

TEST_P(anchor_drop_test, drop_right_short) {
  m_options.min_overlap = tseq("yz").size();
  m_scaffold.add(0, tseq("abcdefghijklmnopqrstuvwxyz"));
  assembly a;
  a.seq = tseq("abcd") + dna_G + tseq("yz");
  a.left_offset = 0;
  a.left_anchor_len = tseq("abcd").size();
  a.right_offset = a.seq.size();
  a.right_anchor_len = 0;
  align(a);

  ASSERT_THAT(m_aligned, ElementsAre(AlignedVarsAre(VarIs(
                             tseq("abcd").size(), dna_G,
                             (tseq("abcdefghijklmnopqrstuvwx").size())))));
  EXPECT_EQ(m_aligned[0].left_offset, 0);
  EXPECT_EQ(m_aligned[0].right_offset,
            tseq("abcdefghijklmnopqrstuvwxyz").size());
  EXPECT_EQ(m_aligned[0].right_anchor_len, tseq("yz").size());
}

INSTANTIATE_TEST_CASE_P(anchor_drop_tests, anchor_drop_test,
                        ::testing::Values(1, 2, 3, 4, 20, 31) /* kmer sizes */);

TEST_F(align_test, two_snps) {
  dna_sequence orig = dna_T + dna_A + dna_T;
  dna_sequence variant = dna_G + dna_A + dna_G;

  m_scaffold.add(5, tseq("abc") + orig + tseq("def"));
  assembly a;
  a.left_offset = 5;
  a.left_anchor_len = 1;
  a.right_offset = 5 + (tseq("abc") + orig + tseq("def")).size();
  a.right_anchor_len = 1;
  a.seq = tseq("abc") + variant + tseq("def");
  align(a);

  EXPECT_THAT(m_aligned, ElementsAre(AlignedVarsAre(
                             VarIs(5 + tseq("abc").size(), dna_G,
                                   5 + (tseq("abc") + dna_T).size()),

                             VarIs(5 + (tseq("abc") + dna_T + dna_A).size(),
                                   dna_G, 5 + (tseq("abc") + orig).size()))));

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(
          RefAssemblyIs(5, 5 + tseq("abc").size()),
          AssemblyIs(5 + tseq("abc").size(), dna_G,
                     5 + (tseq("abc") + dna_T).size()),
          RefAssemblyIs(5 + (tseq("abc") + dna_T).size(),
                        5 + (tseq("abc") + dna_T + dna_A).size()),
          AssemblyIs(5 + (tseq("abc") + dna_T + dna_A).size(), dna_G,
                     5 + (tseq("abc") + orig).size()),
          RefAssemblyIs(5 + (tseq("abc") + orig).size(),
                        5 + (tseq("abc") + orig + tseq("def")).size())));
}

TEST_F(align_test, multi_extent) {
  m_scaffold.add(0, tseq("abc"));
  m_scaffold.add(100, tseq("ghi"));

  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = 1;
  a.right_offset = 100 + tseq("ghi").size();
  a.right_anchor_len = 1;
  a.seq = tseq("abcdefghi");
  align(a);

  EXPECT_THAT(
      m_aligned,
      ElementsAre(AlignedVarsAre(VarIs(tseq("abc").size(), tseq("def"), 100))));

  EXPECT_THAT(m_assemblies,
              ElementsAre(RefAssemblyIs(0, tseq("abc").size()),
                          AssemblyIs(tseq("abc").size(), tseq("def"), 100),
                          RefAssemblyIs(100, 100 + tseq("ghi").size())));
}

TEST_F(align_test, insert_base) {
  m_scaffold.add(0, tseq("abcdef"));

  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = 1;
  a.right_offset = tseq("abcdef").size();
  a.right_anchor_len = 1;
  a.seq = tseq("abc") + dna_T + tseq("def");
  align(a);

  EXPECT_THAT(m_aligned, ElementsAre(AlignedVarsAre(VarIs(
                             tseq("abc").size(), dna_T, tseq("abc").size()))));

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(RefAssemblyIs(0, tseq("abc").size()),
                  AssemblyIs(tseq("abc").size(), dna_T, tseq("abc").size()),
                  RefAssemblyIs(tseq("abc").size(), tseq("abcdef").size())));
}

TEST_F(align_test, delete_base) {
  m_scaffold.add(0, tseq("abc") + dna_T + tseq("def"));

  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = 1;
  a.right_offset = (tseq("abc") + dna_T + tseq("def")).size();
  a.right_anchor_len = 1;
  a.seq = tseq("abcdef");
  align(a);

  EXPECT_THAT(
      m_aligned,
      ElementsAre(AlignedVarsAre(VarIs(tseq("abc").size(), dna_sequence(),
                                       (tseq("abc") + dna_T).size()))));
  EXPECT_THAT(
      m_assemblies,
      ElementsAre(RefAssemblyIs(0, tseq("abc").size()),
                  AssemblyIs(tseq("abc").size(), dna_sequence(),
                             (tseq("abc") + dna_T).size()),
                  RefAssemblyIs((tseq("abc") + dna_T).size(),
                                (tseq("abc") + dna_T + tseq("def")).size())));
}

TEST_F(align_test, sorts_by_left) {
  m_scaffold.add(0, tseq("abc") + dna_T + tseq("def") + dna_T + tseq("ghi") +
                        dna_T + tseq("jkl"));

  std::vector<assembly> asms;

  assembly a;
  a.left_offset = 0;
  a.left_anchor_len = 1;
  a.right_offset =
      (tseq("abc") + dna_T + tseq("def") + dna_T + tseq("ghi")).size();
  a.right_anchor_len = 1;
  a.seq = tseq("abc") + dna_C + tseq("def") + dna_C + tseq("ghi");

  asms.push_back(a);

  a.left_offset = (tseq("abc") + dna_T + tseq("def")).size();
  a.left_anchor_len = 1;
  a.right_offset = (tseq("abc") + dna_T + tseq("def") + dna_T + tseq("ghi") +
                    dna_T + tseq("jk"))
                       .size();
  a.right_anchor_len = 1;
  a.seq = tseq("def") + dna_C + tseq("ghi") + dna_C + tseq("jkl");
  asms.push_back(a);

  expect_sorted(assembly::left_offset_less_than);
}

// Examples that have been seen in the wild, and where they should align.
class wild_align_test : public align_test {
 public:
  std::pair<dna_sequence, dna_sequence> p(const std::string& ref,
                                          const std::string& var) {
    return std::make_pair(dna_sequence(ref), dna_sequence(var));
  };

  void SetUp() override { disable_test_sequence_expansion(); }

  void TearDown() override { enable_test_sequence_expansion(); }

  void wild_align_seq(
      std::vector<std::pair<dna_sequence, dna_sequence>> ref_vs_var) {
    ASSERT_GT(m_pre_ref.size(), 0);
    ASSERT_GT(m_post_ref.size(), 0);

    ASSERT_EQ(m_aligned.size(), 0);

    ref_vs_var.insert(ref_vs_var.begin(), std::make_pair(m_pre_ref, m_pre_ref));
    ref_vs_var.push_back(std::make_pair(m_post_ref, m_post_ref));

    dna_sequence ref_seq;
    dna_sequence var_seq;
    for (const auto& v : ref_vs_var) {
      const auto& ref = v.first;
      const auto& var = v.second;

      if (ref != var) {
        m_expected_vars.push_back(
            VarIs(ref_seq.size(), var, ref_seq.size() + ref.size()));
      }
      ref_seq += ref;
      var_seq += var;
    }

    m_scaffold.add(0, ref_seq);

    assembly a;
    a.seq = var_seq;
    a.left_offset = 0;
    a.left_anchor_len = m_pre_ref.size();
    a.right_offset = ref_seq.size();
    a.right_anchor_len = m_post_ref.size();
    align(a);
  }

  dna_sequence m_pre_ref;
  dna_sequence m_post_ref;
  std::vector<Matcher<const aligned_var&>> m_expected_vars;
};

TEST_F(wild_align_test, aligns_smaller_near_ends) {
  m_pre_ref = dna_sequence{
      "TGTGATAGTAGATCAGAGGTGGAGACATCAACGTAAACTTATGTTTAGTTTAATATAGACACACACAGTTCT"
      "ACATAGAAAACTTTATAATTAGGTGTGT"};
  m_post_ref = dna_sequence{
      "ACCATTCTGAAAGGAATCAGGCTCTTTGAAGAAATGTCTGATACTAGAACTGGGACAGTAAATATAGGAGCC"
      "AGGATAATCTGGAAGTATCAGAAAGTAAGTACTAAAAAAATTAAAATATATCAAACAAAAATAAAAGCCAAT"
      "AAAA"};

  wild_align_seq({p("G", "A"),                                //
                  p("TAGGTAGGTTAGACAC", "TAGGTAGGTTAGACAC"),  //
                  p("G", "A"),
                  p("CACATATACTTCCTAGCATTGC", "CACATATACTTCCTAGCATTGC"),  //
                  p("T", "C"),
                  p("AATGAGGGA",  //
                    "AATGAGGGA"),
                  p("C", "G"),
                  p("AAGATACAATGTGC",  //
                    "AAGATACAATGTGC"),
                  p("", "TC"),
                  p("ATTCAGCAGCCA",  //
                    "ATTCAGCAGCCA"),
                  p("C", "G"),
                  p("ATGTAAGTTTTCC",  //
                    "ATGTAAGTTTTCC"),
                  p("C", "T")});
  EXPECT_THAT(m_aligned, ElementsAre(AlignedVarsAreArray(m_expected_vars)));
}

}  // namespace variants

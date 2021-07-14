#include "modules/variants/calc_coverage.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;
using namespace coverage_testutil;

class calc_coverage_test : public assemble_test {
 public:
  calc_coverage_test() { m_options.penalize_directional_coverage = false; }
  void start_calc() { m_calc.emplace(m_options, test_output()); }
  void add(assembly a) { m_calc->add(make_unique<assembly>(a)); }
  void flush() {
    m_calc.reset();
    expect_sorted(assembly::left_offset_less_than);
  }

 protected:
  boost::optional<calc_coverage> m_calc;
};

TEST_F(calc_coverage_test, single) {
  use_reads({tseq("bcdefg")});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("abcdefgh");
  a.right_offset = 1;
  add(a);
  flush();

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(AllOf(CoverageIs(0 + over("a", 0) + 0 + over("bcdefg", 1) + 0 + over("h", 0) + 0),
                        PairCoverageIs(0 + over("abcdefgh", 0) + 0))));
}

TEST_F(calc_coverage_test, dual) {
  use_reads({tseq("bcdefg"), tseq("defgh")});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("abcdefghi");
  a.right_offset = 1;
  a.right_pair_matches.push_back(read_id_for_seq(tseq_rc("defgh")));
  add(a);
  flush();

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(
          AllOf(CoverageIs(0 + over("a", 0) + 0 + over("bc", 1) + 1 + over("defg", 2) + 1 +
                           over("h", 1) + 0 + over("i", 0) + 0),
                PairCoverageIs(0 + over("abc", 0) + 0 + over("defgh", 1) + 0 + over("i", 0) + 0))));
}

TEST_F(calc_coverage_test, left_edge) {
  use_reads({tseq("abc")});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("abcdefg");
  a.right_offset = 1;

  add(a);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(CoverageIs(0 + over("abc", 1) + 0 + over("defg", 0) + 0)));
}

TEST_F(calc_coverage_test, right_edge) {
  use_reads({tseq("efg")});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("abcdefg");
  a.right_offset = 1;

  add(a);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(CoverageIs(0 + over("abcd", 0) + 0 + over("efg", 1) + 0)));
}

TEST_F(calc_coverage_test, exact_full) {
  use_reads({tseq("abcdefg")});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("abcdefg");
  a.seq = tseq("abcdefg");
  a.right_offset = 1;
  a.left_pair_matches.push_back(read_id_for_seq(tseq_rc("abcdefg")));

  add(a);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(CoverageIs(0 + over("abcdefg", 1) + 0),
                                              PairCoverageIs(0 + over("abcdefg", 1) + 0))));
}

TEST_F(calc_coverage_test, split) {
  use_reads({tseq("abcdef"), tseq("cdefg")});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("abcdefg");
  a.right_offset = 1000;
  add(a);
  flush();

  ASSERT_THAT(m_assemblies, ElementsAre(CoverageIs(0 + over("ab", 1) + 1 + over("cdef", 2) + 1 +
                                                   over("g", 1) + 0)));

  auto split_pair = split_assembly(make_unique<assembly>(m_assemblies[0]), tseq("ab").size(), 500);
  EXPECT_THAT(split_pair, Pair(Pointee(CoverageIs(0 + over("ab", 1) + 1)),
                               Pointee(CoverageIs(1 + over("cdef", 2) + 1 + over("g", 1) + 0))));

  split_pair = split_assembly(make_unique<assembly>(m_assemblies[0]), tseq("ab").size() + 1, 500);
  EXPECT_THAT(split_pair, Pair(Pointee(CoverageIs(0 + over("ab", 1) + 1 + 2)),
                               Pointee(CoverageIs(over("cdef", 2) + 1 + over("g", 1) + 0))));
}

TEST_F(calc_coverage_test, long_overlap) {
  dna_sequence seq = tseq("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
  dna_sequence reads_seq = tseq("XXXX") + seq + tseq("XXXX");
  use_ref_parts({{0, reads_seq}});
  std::vector<dna_sequence> reads;
  for (unsigned i = 1; i + tseq("XXXX").size() + 1 <= reads_seq.size(); i += tseq("X").size()) {
    reads.emplace_back(reads_seq.subseq(i, tseq("XXXX").size() + 1));
  }

  use_reads(reads);

  start_calc();
  assembly a;
  a.left_offset = tseq("XXXX").size();
  a.seq = seq;
  a.right_offset = reads_seq.size() - tseq("XXXX").size();
  add(a);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(CoverageIs(rpt(seq.size() + 1, 4))));
}

Matcher<assembly> OtherDepthIs(unsigned expected_other_depth) {
  return Field(&assembly::other_depth, expected_other_depth);
}

Matcher<assembly> RefDepthIs(unsigned expected_ref_depth) {
  return Field(&assembly::ref_depth, expected_ref_depth);
}

TEST_F(calc_coverage_test, basic_with_ref) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd") + dna_A + tseq("efg")});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcd").size();
  a.seq = dna_A;
  a.right_offset = tseq("abcd").size() + 1;
  add(a);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(CoverageIs(cov(1) + 1), OtherDepthIs(0))));
}

TEST_F(calc_coverage_test, basic_with_ref_left) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd") + dna_A, tseq("cd") + dna_A + tseq("efg")});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcd").size();
  a.seq = dna_A;
  a.right_offset = tseq("abcd").size() + 1;
  add(a);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(CoverageIs(cov(2) + 1), OtherDepthIs(0))));
}

TEST_F(calc_coverage_test, basic_with_ref_right) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd") + dna_A + tseq("efg"), dna_A + tseq("efgh")});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcd").size();
  a.seq = dna_A;
  a.right_offset = tseq("abcd").size() + 1;
  add(a);
  flush();

  EXPECT_THAT(m_assemblies,
              ElementsAre(AllOf(CoverageIs(cov(1) + 2), OtherDepthIs(0), RefDepthIs(0))));
}

TEST_F(calc_coverage_test, basic_with_ref_coverage) {
  use_ref_parts({{0, tseq("abcd") + dna_T + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd") + dna_T + tseq("efg")});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcd").size();
  a.seq = dna_A;
  a.right_offset = tseq("abcd").size() + 1;
  add(a);
  flush();

  EXPECT_THAT(m_assemblies,
              ElementsAre(AllOf(CoverageIs(cov(0) + 0), OtherDepthIs(0), RefDepthIs(1))));
}

TEST_F(calc_coverage_test, insert_between_scaffolds) {
  use_ref_parts({{0, tseq("abcd")}, {1000, tseq("ijkl")}});
  use_reads({tseq("abcd"), tseq("cdef"), tseq("efgh"), tseq("ghij"), tseq("ijkl")});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcd").size();
  a.seq = tseq("efgh");
  a.right_offset = 1000;
  add(a);
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(AllOf(CoverageIs(1 + over("ef", 2) + 1 + over("gh", 2) + 1),
                                              OtherDepthIs(0), RefDepthIs(0))));
}

class coverage_extent_edge_test
    : public calc_coverage_test,
      public WithParamInterface<std::tuple<
          int /* edge offset */,
          bool /* true if we add a second assembly to "skip" past the end of the extent instead of just closing out the calculator */>> {
 public:
  coverage_extent_edge_test()
      : m_edge_offset(std::get<0>(GetParam())), m_second_assembly(std::get<1>(GetParam())) {}

 protected:
  const int m_edge_offset;
  const bool m_second_assembly;
};

constexpr int k_extent_edge_read_len = 40;

// Tests edge cases for when the region of reference we need to
// calculate coverage for overlaps the trailing edge of an extent.
TEST_P(coverage_extent_edge_test, part_overlaps_edge) {
  dna_sequence next_ref = tseq("ghijklmnopqrstuvwxyz");
  ASSERT_LT(m_edge_offset, next_ref.size());
  dna_sequence main_ref = tseq("abcdef") + next_ref.subseq(0, m_edge_offset);
  use_ref_parts({{0, main_ref}, {1000, tseq("ijkl")}});
  int last_ref_read_pos = main_ref.size() - k_extent_edge_read_len;
  dna_sequence last_ref_read = main_ref.subseq(last_ref_read_pos, k_extent_edge_read_len);
  use_reads({tseq("abcd"), tseq("cdef"), tseq("efgh"), last_ref_read});

  ASSERT_EQ(k_extent_edge_read_len, m_options.seqset->max_read_len());

  start_calc();
  assembly a;
  a.seq = dna_T;
  int var_start_pos = tseq("abcdef").size();
  a.left_offset = var_start_pos;
  a.right_offset = var_start_pos + 1;
  add(a);

  if (m_second_assembly) {
    // Add another assembly in the next extent so we have to skip to the
    // next extent instead of just closing out the calculation.
    a.left_offset = 1000 + tseq("ij").size();
    a.right_offset = 1000 + tseq("ij").size() + 1;
    add(a);
  }
  flush();

  int expected_ref_depth = 1;
  if (last_ref_read_pos < var_start_pos) {
    // last_ref_read includes ref coverage for the first variant.
    ++expected_ref_depth;
  }

  auto first_assembly_matcher =
      AllOf(CoverageIs(cov(0) + 0), OtherDepthIs(0), RefDepthIs(expected_ref_depth));
  if (m_second_assembly) {
    auto second_assembly_matcher = AllOf(CoverageIs(cov(0) + 0), OtherDepthIs(0), RefDepthIs(0));
    EXPECT_THAT(m_assemblies, ElementsAre(first_assembly_matcher, second_assembly_matcher));
  } else {
    EXPECT_THAT(m_assemblies, ElementsAre(first_assembly_matcher));
  }
}

INSTANTIATE_TEST_CASE_P(coverage_extent_edge_tests, coverage_extent_edge_test,
                        ::testing::Combine(
                            // Test for a bunch of positions right around the edge of the extent
                            ::testing::Range(k_extent_edge_read_len - 2,
                                             k_extent_edge_read_len + 4),
                            // Test for both with and without an assembly in a second extent, so
                            // the result is the same whether we skip to the next extent or whether
                            // we just close down the calculator.
                            ::testing::Bool()));

class max_coverage_paths_test : public calc_coverage_test,
                                public WithParamInterface<int /* max coverage paths */> {
 public:
  max_coverage_paths_test() { m_options.max_coverage_paths = GetParam(); }
};

TEST_P(max_coverage_paths_test, snp_with_other_coverage) {
  use_ref_parts({{0, tseq("abcd") + dna_G + tseq("efghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd") + dna_T + tseq("efg"), tseq("abcd") + dna_A + tseq("efg"),
             tseq("bcd") + dna_A + tseq("efgh")});

  start_calc();
  assembly a1;
  a1.assembly_id = 1;
  a1.left_offset = tseq("abcd").size();
  a1.seq = dna_T;
  a1.right_offset = tseq("abcd").size() + 1;
  add(a1);

  assembly a2;
  a2.assembly_id = 2;
  a2.left_offset = tseq("abcd").size();
  a2.seq = dna_A;
  a2.right_offset = tseq("abcd").size() + 1;
  add(a2);
  flush();

  constexpr int k_broken_other_depth = 0;  // TODO(nils): This should be 1?

  auto expected = UnorderedElementsAre(
      AllOf(AssemblyIdIs(1), CoverageIs(rpt(2, 1)), OtherDepthIs(2), RefDepthIs(0)),
      AllOf(AssemblyIdIs(2), CoverageIs(rpt(2, 2)), OtherDepthIs(k_broken_other_depth),
            RefDepthIs(0)));
  if (m_options.max_coverage_paths < 2) {
    EXPECT_THAT(m_assemblies, Not(expected));
  } else {
    EXPECT_THAT(m_assemblies, expected);
  }
}

INSTANTIATE_TEST_CASE_P(max_coverage_paths_succeeds, max_coverage_paths_test, ::testing::Values(3));
INSTANTIATE_TEST_CASE_P(max_coverage_paths_skips_ref, max_coverage_paths_test,
                        ::testing::Values(2));
INSTANTIATE_TEST_CASE_P(max_coverage_paths_fails, max_coverage_paths_test, ::testing::Values(1));

class calc_coverage_stress_test : public calc_coverage_test, public WithParamInterface<int> {
 public:
  calc_coverage_stress_test() { m_num_vars = GetParam(); }

 protected:
  int m_num_vars = 0;
  static constexpr int k_max_count = 5;
};

constexpr int calc_coverage_stress_test::k_max_count;

TEST_P(calc_coverage_stress_test, non_compound_heterozygous) {
  static const std::string k_alpha = "abcdefghijklmnopqrstuvwxyz";
  CHECK_LT(m_num_vars, k_alpha.size());

  m_options.min_anchor_drop_overlap = tseq("a").size();

  size_t seed = std::random_device()();
  std::cout << "Using random seed " << seed << "\n";
  std::mt19937 rand_source(seed);

  std::vector<bool> is_delete(m_num_vars, false);
  std::vector<bool> is_insert(m_num_vars, false);
  std::vector<int> var_count(m_num_vars, 0);
  std::vector<int> pair_count(m_num_vars, 0);
  std::vector<int> ref_count(m_num_vars, 0);
  std::vector<assembly> variants;
  int total_count = 0;
  std::uniform_int_distribution<int> rand_bool(0, 1);
  std::uniform_int_distribution<int> rand_count(1, k_max_count);
  std::uniform_int_distribution<int> rand_offset(0, tseq("a").size());
  for (int i = 0; i < m_num_vars; ++i) {
    is_delete[i] = rand_bool(rand_source);
    if (is_delete[i]) {
      is_insert[i] = false;
    } else {
      is_insert[i] = rand_bool(rand_source);
    }
    var_count[i] = rand_count(rand_source);
    total_count += var_count[i];
  }

  dna_sequence ref_seq = tseq("X");
  dna_sequence all_var_seq = tseq("X");
  for (int i = 0; i < m_num_vars; ++i) {
    assembly a;
    a.assembly_id = i + 1000;
    ref_seq += tseq(k_alpha.substr(i, 1));
    all_var_seq += tseq(k_alpha.substr(i, 1));
    a.left_offset = ref_seq.size();
    if (!is_insert[i]) {
      ref_seq += dna_T;
    }
    if (!is_delete[i]) {
      a.seq += dna_A;
      all_var_seq += dna_A;
    }
    a.right_offset = ref_seq.size();
    variants.push_back(a);
  }
  ref_seq += tseq(k_alpha.substr(m_num_vars, 1));
  ref_seq += tseq("X");
  all_var_seq += tseq(k_alpha.substr(m_num_vars, 1));
  all_var_seq += tseq("X");

  std::cout << "Ref seq:     " << ref_seq << "\n";
  std::cout << "All var seq: " << all_var_seq << "\n";

  use_ref_parts({{0, ref_seq}});

  std::vector<int> var_count_left = var_count;
  std::vector<dna_sequence> reads;
  std::vector<std::set<dna_sequence>> variants_paired_reads;
  variants_paired_reads.resize(m_num_vars);
  while (total_count) {
    std::vector<bool> var_present(m_num_vars, false);
    bool pair_matches = rand_bool(rand_source);

    std::vector<int> var_has_pair;

    for (int i = 0; i < m_num_vars; ++i) {
      if (var_count_left[i] && rand_bool(rand_source)) {
        var_present[i] = true;
        var_count_left[i]--;
        CHECK_GT(total_count, 0);
        --total_count;
      } else {
        var_present[i] = false;
        ++ref_count[i];
      }
    }

    dna_sequence seq = tseq("X");
    for (int i = 0; i < m_num_vars; ++i) {
      seq += tseq(k_alpha.substr(i, 1));
      if (var_present[i]) {
        seq += variants[i].seq;
      } else {
        const auto& a = variants[i];
        seq += ref_seq.subseq(a.left_offset, a.right_offset - a.left_offset);
      }
    }
    seq += tseq(k_alpha.substr(m_num_vars, 1));
    seq += tseq("X");
    // Make this read be not perfectly aligned, so we have a bunch of different reads.
    // TODO(nils): Trim a random offset off left and right.  But we'll have to get
    // variable length read coverage working right first...
    aoffset_t offset = rand_offset(rand_source);
    aoffset_t len = seq.size() - tseq("a").size();
    CHECK_GT(len, 0);
    CHECK_LE(offset + len, seq.size());
    seq = seq.subseq(offset, len);
    reads.push_back(seq);

    if (pair_matches) {
      for (int i = 0; i < m_num_vars; ++i) {
        if (var_present[i]) {
          if (!variants_paired_reads[i].count(seq)) {
            variants_paired_reads[i].insert(seq);
            pair_count[i]++;
          }
        }
      }
    }
    if (reads.size() < 15) {
      std::cout << "Adding read " << seq;
      if (pair_matches) {
        std::cout << ", PAIR suport\n";
      } else {
        std::cout << ", no pair support\n";
      }
    }
  }
  std::cout << "Generated " << reads.size() << " reads for tree test\n";
  use_reads(reads);

  m_options.max_coverage_paths = 1ULL << std::min(m_num_vars, 60);
  start_calc();
  CHECK_EQ(variants.size(), m_num_vars);
  for (int i = 0; i < m_num_vars; ++i) {
    assembly a = variants[i];
    for (dna_sequence seq : variants_paired_reads[i]) {
      uint32_t read_id = read_id_for_seq(seq.rev_comp());
      if (rand_bool(rand_source)) {
        a.left_pair_matches.push_back(read_id);
      } else {
        a.right_pair_matches.push_back(read_id);
      }
    }
    std::cout << "Adding variant: " << a << "\n";
    check_assembly(a, "calc_coverage_test");
    add(a);
  }
  flush();

  std::vector<Matcher<assembly>> matchers;
  for (int i = 0; i < m_num_vars; ++i) {
    matchers.push_back(
        AllOf(Field(&assembly::assembly_id, i + 1000),
              CoverageIs(is_delete[i] ? rpt(1, var_count[i]) : rpt(2, var_count[i])),
              PairCoverageIs(is_delete[i] ? rpt(1, pair_count[i]) : rpt(2, pair_count[i])),
              RefDepthIs(ref_count[i])));
  }
  EXPECT_THAT(m_assemblies, ElementsAreArray(matchers));
}

INSTANTIATE_TEST_CASE_P(calc_coverage_stress_tests, calc_coverage_stress_test,
                        ::testing::Values(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20));

}  // namespace

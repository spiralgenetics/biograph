#include "modules/variants/read_cov.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;
using namespace coverage_testutil;

Matcher<const std::vector<read_coverage_read_t>&> ReadCoverageIs(
    std::vector<read_coverage_read_t> cov) {
  return ContainerEq(cov);
}

const int deref_cov_asm_len(const assembly& a) { return a.read_coverage->assembly_len(); }

const std::vector<read_coverage_read_t>& deref_cov(const assembly& a) {
  return a.read_coverage->reads();
}

class read_cov_test : public assemble_test, public WithParamInterface<bool> {
 public:
  Matcher<const assembly&> AsmReadCoverageIs(
      const dna_sequence& seq, std::set<std::pair<dna_sequence, aoffset_t>> expected_seqs) {
    read_coverage_set expected;
    for (const auto& expected_seq : expected_seqs) {
      const auto& seq = expected_seq.first;
      aoffset_t offset = expected_seq.second;

      for (uint32_t read_id : get_read_ids(seq)) {
        expected.insert(offset, read_id, m_options.readmap->get_readlength(read_id));
      }
    }

    read_coverage_t expected_cov = expected.build_and_clear(seq.size());
    return AllOf(Field(&assembly::seq, Eq(seq)), ResultOf(deref_cov_asm_len, Eq(seq.size())),
                 ResultOf(deref_cov, ReadCoverageIs(expected_cov.reads())));
  }

  read_cov_test() { m_rev_comp = GetParam(); }

  void use_reference(const std::string& reference_dir, const std::string& scaffold_name) {
    assemble_test::use_reference(reference_dir, scaffold_name);

    m_ref_end_pos = m_options.scaffold->end_pos();

    if (m_rev_comp) {
      m_scaffold = m_scaffold.rev_comp();
      m_options.scaffold = &m_scaffold;
    }
  }

  void use_read_ref(std::vector<std::pair<aoffset_t, dna_sequence>> parts) {
    use_ref_parts(parts);

    m_ref_end_pos = m_options.scaffold->end_pos();

    if (m_rev_comp) {
      m_scaffold = m_scaffold.rev_comp();
      m_options.scaffold = &m_scaffold;
    }
  }
  void start_calc() {
    auto cov = make_unique<read_cov>(m_options, test_output());
    m_ref_adder.emplace(m_options, m_options.seqset->max_read_len() * 2, false /* whole ref */,
                        0 /* length limit */, std::move(cov));
  }
  void add(assembly a) {
    if (m_rev_comp) {
      rev_asm(a);
    }
    m_ref_adder->add(make_unique<assembly>(a));
  }
  std::unordered_set<uint32_t> get_read_ids(const dna_sequence& seq) {
    seqset_range r = m_options.seqset->find(seq);
    CHECK(r.valid()) << seq;
    std::unordered_set<uint32_t> read_ids;
    for (const auto& read : m_options.readmap->get_prefix_reads(r)) {
      if (read.size() != int(seq.size())) {
        continue;
      }
      read_ids.insert(read.get_read_id());
    }
    return read_ids;
  }
  uint32_t get_read_id(const dna_sequence& seq) {
    auto ids = get_read_ids(seq);
    CHECK_EQ(ids.size(), 1) << seq;
    return *ids.begin();
  }
  void pad(std::ostream& os, int n) {
    while (n > 0) {
      os << " ";
      --n;
    }
  }
  template <typename C>
  std::string print_all(const C& c) {
    std::stringstream os;
    for (const assembly& a : c) {
      os << "Assembly: " << a << "\nRead coverage: " << print_read_cov(a) << "\n";
    }
    return os.str();
  }
  std::string print_read_cov(const assembly& a) {
    if (!a.read_coverage) {
      return "(no coverage)";
    }
    std::stringstream os;

    const read_coverage_t& cov = *a.read_coverage;

    if (cov.reads().empty()) {
      return "(empty coverage)";
    }

    aoffset_t npad = 0 - cov.reads().begin()->offset;
    if (npad < 0) {
      npad = 0;
    }
    os << "Assembly with " << cov.reads().size() << " reads:\n";
    pad(os, npad);
    os << "\"" << a.seq.as_string() << "\"\n";
    for (const auto& rd : cov.reads()) {
      for (uint32_t read_id : rd.read_ids) {
        pad(os, npad + rd.offset);
        os << m_options.readmap->get_read_by_id(read_id).get_seqset_entry().sequence() << " (@"
           << rd.offset << ")\n";
      }
    }

    return os.str();
  }

  void flush_and_check() {
    m_ref_adder.reset();
    expect_sorted(assembly::left_offset_less_than);

    if (m_rev_comp) {
      for (auto* collection : {&m_assemblies, &m_ref_assemblies, &m_non_ref_assemblies}) {
        for (auto& a : *collection) {
          rev_asm(a);
        }
        std::reverse(collection->begin(), collection->end());
      }
      m_scaffold = m_scaffold.rev_comp();
    }

    for (auto* collection : {&m_assemblies, &m_ref_assemblies, &m_non_ref_assemblies}) {
      for (auto& a : *collection) {
        verify_assembly(a);
      }
    }
  }

  void verify_assembly(const assembly& a) {
    if (!a.read_coverage) {
      return;
    }

    const read_coverage_t& cov = *a.read_coverage;
    EXPECT_EQ(cov.assembly_len(), a.seq.size());

    for (const auto& cov_entry : cov.reads()) {
      for (uint32_t read_id : cov_entry.read_ids) {
        readmap::read rd = m_options.readmap->get_read_by_id(read_id);
        aoffset_t read_start_offset = cov_entry.offset;
        aoffset_t read_end_offset = read_start_offset + rd.size();

        EXPECT_GT(read_end_offset, 0) << a << "\ncoverage:\n" << print_read_cov(a);
        EXPECT_LT(read_start_offset, aoffset_t(a.seq.size())) << a << "\ncoverage:\n"
                                                              << print_read_cov(a);
        if (read_end_offset <= 0 || read_start_offset >= aoffset_t(a.seq.size())) {
          return;
        }

        aoffset_t trunc_start = std::max<aoffset_t>(read_start_offset, 0);
        aoffset_t trunc_end = std::min<aoffset_t>(read_end_offset, a.seq.size());
        aoffset_t trunc_len = trunc_end - trunc_start;
        EXPECT_GT(trunc_len, 0);
        if (trunc_len <= 0) {
          break;
        }
        dna_sequence seq_covered = a.seq.subseq(trunc_start, trunc_len);
        dna_sequence read_seq = rd.get_seqset_entry().sequence();
        dna_sequence read_seq_covering =
            read_seq.subseq(trunc_start - read_start_offset, trunc_len);

        EXPECT_EQ(seq_covered, read_seq_covering) << a << "\ncoverage:\n" << print_read_cov(a);
      }
    }

    auto overlaps = cov.get_overlaps();
    auto overlap_min_max = cov.get_overlap_min_max();
    if (overlaps.empty()) {
      EXPECT_THAT(overlap_min_max, Pair(0, 0));
    } else {
      EXPECT_THAT(overlap_min_max, Pair(*std::min_element(overlaps.begin(), overlaps.end()),
                                        *std::max_element(overlaps.begin(), overlaps.end())));
    }
  }

  void rev_asm(assembly& a) {
    reverse_assembly_in_place(&a, m_options.readmap, m_options.scaffold->end_pos());
  }

 protected:
  bool m_rev_comp;
  scaffold m_rev_scaffold;
  boost::optional<add_ref> m_ref_adder;
  aoffset_t m_ref_end_pos;
};

TEST_P(read_cov_test, var_cov) {
  use_read_ref({{0, tseq("abcdefghijklmnopqrs")}});
  use_reads({tseq("defgh"), tseq("efgh") + dna_A, tseq("IJKL"), dna_A + tseq("IJKL") + dna_A,
             dna_A + tseq("mnop"), tseq("mnopq")});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcdefgh").size();
  a.seq = dna_A + tseq("IJKL") + dna_A;
  a.right_offset = tseq("abcdefghijkl").size();
  add(a);
  flush_and_check();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AsmReadCoverageIs(dna_A + tseq("IJKL") + dna_A,
                                            {{tseq("IJKL"), 1},
                                             {dna_A + tseq("IJKL") + dna_A, 0},
                                             {tseq("efgh") + dna_A, -tseq("efgh").size()},
                                             {dna_A + tseq("mnop"), 1 + tseq("IJKL").size()}})))
      << print_all(m_non_ref_assemblies);

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(                             //
                  AsmReadCoverageIs(tseq("abcdefgh"),  //
                                    {{tseq("defgh"), tseq("abc").size()},
                                     {tseq("efgh") + dna_A, tseq("abcd").size()}}),
                  AsmReadCoverageIs(tseq("ijkl"),  //
                                    {}),
                  AsmReadCoverageIs(tseq("mnopqrs"),  //
                                    {{dna_A + tseq("mnop"), -1}, {tseq("mnopq"), 0}})))
      << print_all(m_ref_assemblies);
}

TEST_P(read_cov_test, ref_cov) {
  use_read_ref({{0, tseq("abcdefgh") + dna_A + tseq("IJKL") + dna_A + tseq("mnopqrs")}});
  use_reads({tseq("defgh"), tseq("efgh") + dna_A, tseq("IJKL"), dna_A + tseq("IJKL") + dna_A,
             dna_A + tseq("mnop"), tseq("mnopq")});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcdefgh").size();
  a.seq = tseq("ijkl");
  a.right_offset = tseq("abcdefgh").size() + 1 + tseq("IJKL").size() + 1;
  add(a);
  flush_and_check();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmReadCoverageIs(tseq("ijkl"), {})))
      << print_all(m_non_ref_assemblies);

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(                             //
                  AsmReadCoverageIs(tseq("abcdefgh"),  //
                                    {{tseq("defgh"), tseq("abc").size()},
                                     {tseq("efgh") + dna_A, tseq("abcd").size()}}),
                  AsmReadCoverageIs(dna_A + tseq("IJKL") + dna_A,  //
                                    {{tseq("efgh") + dna_A, -tseq("efgh").size()},
                                     {dna_A + tseq("IJKL") + dna_A, 0},
                                     {tseq("IJKL"), 1},
                                     {dna_A + tseq("mnop"), 1 + tseq("IJKL").size()}}),
                  AsmReadCoverageIs(tseq("mnopqrs"),  //
                                    {{dna_A + tseq("mnop"), -1}, {tseq("mnopq"), 0}})))
      << print_all(m_ref_assemblies);
}

TEST_P(read_cov_test, inserts) {
  use_read_ref({{0, tseq("abcdefghijklmnopqrs")}});
  use_reads({
      // Coverage for first insert:
      tseq("efg") + dna_T + tseq("hij"),
      // Coverage for second insert:
      tseq("defg") + dna_A + tseq("hi"),
      // Coverage for both inserts in the same allele; these should not be included.
      tseq("efg") + dna_A + dna_T + tseq("hij"),
      tseq("efg") + dna_T + dna_A + tseq("hij"),
  });

  start_calc();
  assembly a;
  a.left_offset = tseq("abcdefg").size();
  a.seq = dna_T;
  a.right_offset = a.left_offset;
  add(a);
  a.seq = dna_A;
  add(a);
  flush_and_check();

  EXPECT_THAT(
      m_non_ref_assemblies,
      UnorderedElementsAre(
          AsmReadCoverageIs(dna_T, {{tseq("efg") + dna_T + tseq("hij"), -tseq("efg").size()}}),
          AsmReadCoverageIs(dna_A, {{tseq("defg") + dna_A + tseq("hi"), -tseq("defg").size()}})))
      << print_all(m_non_ref_assemblies);
}

namespace {

std::vector<dna_sequence> make_reads_multiplied(
    std::vector<std::pair<int, dna_sequence>> reads_and_counts) {
  std::vector<dna_sequence> result;

  for (const auto& r_and_c : reads_and_counts) {
    for (int i = 0; i < r_and_c.first; ++i) {
      result.push_back(r_and_c.second);
    }
  }
  return result;
}

}  // namespace

TEST_P(read_cov_test, calc_depths) {
  use_read_ref({{0, tseq("abcd") + dna_A + dna_T + dna_A + tseq("ijklmnopqrs")}});

  use_reads(make_reads_multiplied({
      {1, dna_A + dna_C + dna_A + tseq("ijkl")},
      {2, (dna_C + dna_A + tseq("ijklm")).rev_comp()},
      {4, dna_A + tseq("ijkl")},

      {8, tseq("bcd") + dna_A},
      {16, tseq("bcd") + dna_A + dna_C},
      {32, tseq("bcd") + dna_A + dna_C + dna_A},
  }));

  start_calc();
  assembly a;
  a.left_offset = (tseq("abcd") + dna_A).size();
  a.seq = dna_C;
  a.right_offset = (tseq("abcd") + dna_A + dna_T).size();
  add(a);
  flush_and_check();

  ASSERT_THAT(m_non_ref_assemblies, SizeIs(1));
  const auto& result_a = m_non_ref_assemblies.front();

  ASSERT_TRUE(result_a.read_coverage);
  const auto& cov = *result_a.read_coverage;

  EXPECT_THAT(cov.calc_depths(true /* include fwd */, true /* include rev */, true /* interbase */,
                              m_options.readmap),
              ElementsAre(1 + 16 + 32, 1 + 2 + 32))
      << print_all(m_non_ref_assemblies);
  EXPECT_THAT(cov.calc_depths(false /* include fwd */, true /* include rev */, true /* interbase */,
                              m_options.readmap),
              ElementsAre(0 + 0 + 0, 0 + 2 + 0))
      << print_all(m_non_ref_assemblies);
  EXPECT_THAT(cov.calc_depths(true /* include fwd */, false /* include rev */, true /* interbase */,
                              m_options.readmap),
              ElementsAre(1 + 16 + 32, 1 + 0 + 32))
      << print_all(m_non_ref_assemblies);

  EXPECT_THAT(cov.calc_depths(true /* include fwd */, true /* include rev */, false /* interbase */,
                              m_options.readmap),
              ElementsAre(1 + 2 + 16 + 32))
      << print_all(m_non_ref_assemblies);
  EXPECT_THAT(cov.calc_depths(false /* include fwd */, true /* include rev */,
                              false /* interbase */, m_options.readmap),
              ElementsAre(0 + 2 + 0 + 0))
      << print_all(m_non_ref_assemblies);
  EXPECT_THAT(cov.calc_depths(true /* include fwd */, false /* include rev */,
                              false /* interbase */, m_options.readmap),
              ElementsAre(1 + 0 + 16 + 32))
      << print_all(m_non_ref_assemblies);

  EXPECT_EQ(cov.get_max_flank(0), 2);
  EXPECT_EQ(cov.get_max_flank(1), 2);

  auto filtered_cov = cov.get_reads_spanning_offset(1 + tseq("ijkl").size());
  EXPECT_THAT(filtered_cov.calc_depths(), ElementsAre(1, 1 + 2));
  EXPECT_THAT(filtered_cov.get_overlaps(),
              ElementsAre(  //
                  (dna_C + dna_A + tseq("ijkl")).size(),
                  // Multiple reads on this one, so it has a overlap of the full read length:
                  (dna_C + dna_A + tseq("ijklm")).size()));
  EXPECT_EQ(filtered_cov.get_max_flank(0), 1);
  EXPECT_EQ(filtered_cov.get_max_flank(1), 2);
}

TEST_P(read_cov_test, wild_inserts) {
  use_biograph("datasets/lambdaToyData/benchmark/proband_lambda.bg");
  use_reference("datasets/lambdaToyData/benchmark/ref_lambda", "lambda");

  dna_sequence partial(
      "AAGAACGTTATAGAGAACCTATCTTTCGGGGATGGGCCTATTGCGTCTAACATAGACACTTTAAGGCTAATGAAGTTTGTAGCTAAGACCGC"
      "TGGGGAGTGAATAGCGGGACACGAATGGTCGGGAAGCAAAACGAAACGGAGGATTCTC");
  dna_sequence full(
      "GCCTAGGCGGGAACGTGGGCCATGGTGGCTGCCGCATGTACTGGCGATTGATCCTCCTGCAACCTGAAGGGACGGCCGCGGGAACGTCTCC"
      "GATAATGAAGGCTTGCACTCATATACTATCCAAGCCACGGGTGATACACCCGTGGCACTAAGAACGTTATAGAGAACCTATCTTTCGGGGAT"
      "GGGCCTATTGCGTCTAACATAGACACTTTAAGGCTAATGAAGTTTGTAGCTAAGACCGCTGGGGAGTGAATAGCGGGACACGAATGGTCGGG"
      "AAGCAAAACGAAACGGAGGATTCTC");

  start_calc();
  assembly a;
  a.left_offset = 30278;
  a.seq = full;
  a.right_offset = a.left_offset;
  add(a);
  a.seq = partial;
  add(a);
  flush_and_check();

  std::vector<assembly> actual_asms;

  for (const auto& out_a : m_non_ref_assemblies) {
    if (out_a.seq == full) {
      actual_asms.push_back(out_a);
    }
  }

  ASSERT_THAT(actual_asms, SizeIs(1));
  const auto& actual = actual_asms[0];

  ASSERT_TRUE(actual.read_coverage);

  EXPECT_EQ(actual.read_coverage->get_tot_read_count(), 431) << print_all(m_non_ref_assemblies);
}

INSTANTIATE_TEST_CASE_P(read_cov_fwd_tests, read_cov_test,
                        ::testing::Values(false /* not rev comp */));
INSTANTIATE_TEST_CASE_P(read_cov_rev_tests, read_cov_test, ::testing::Values(true /* rev comp */));

}  // namespace variants

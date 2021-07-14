#include "modules/variants/pair_edge_cov.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/pair_cov.h"
#include "modules/variants/read_cov.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;
using namespace coverage_testutil;

MATCHER_P(ReadIdSetSizeIs, size, "") { return int(arg.size()) == int(size); }

Matcher<const edge_coverage_t&> EdgeCoverageIs(int var_start, int var_end, int ref_start,
                                               int ref_end) {
  return AllOf(Field(&edge_coverage_t::variant_start, ReadIdSetSizeIs(var_start)),
               Field(&edge_coverage_t::variant_end, ReadIdSetSizeIs(var_end)),
               Field(&edge_coverage_t::reference_start, ReadIdSetSizeIs(ref_start)),
               Field(&edge_coverage_t::reference_end, ReadIdSetSizeIs(ref_end)));
}

const edge_coverage_t& deref_edge(const assembly& a) { return *a.edge_coverage; }

Matcher<const assembly&> AsmEdgeCoverIs(int var_start, int var_end, int ref_start, int ref_end) {
  return ResultOf(deref_edge, EdgeCoverageIs(var_start, var_end, ref_start, ref_end));
}

Matcher<const assembly&> AsmEdgeCoverAndIdIs(int var_start, int var_end, int ref_start, int ref_end,
                                             int assembly_id) {
  return AllOf(ResultOf(deref_edge, EdgeCoverageIs(var_start, var_end, ref_start, ref_end)),
               AssemblyIdIs(assembly_id));
}

class pair_edge_cov_test : public assemble_test, public WithParamInterface<bool> {
 public:
  pair_edge_cov_test() {
    m_options.min_pair_distance = 1;
    m_rev_comp = GetParam();
  }

  void use_pair_edge_ref(std::vector<std::pair<aoffset_t, dna_sequence>> parts) {
    use_ref_parts(parts);

    m_ref_end_pos = m_options.scaffold->end_pos();

    if (m_rev_comp) {
      m_scaffold = m_scaffold.rev_comp();
      m_options.scaffold = &m_scaffold;
    }
  }
  void start_calc() {
    auto edge_step = make_unique<pair_edge_cov>(m_options, test_output());
    auto pair_step = make_unique<pair_cov>(m_options, std::move(edge_step));
    auto read_step = make_unique<read_cov>(m_options, std::move(pair_step));
    m_cov.emplace(m_options, m_options.max_pair_distance, false /* whole ref */,
                  0 /* length limit */, std::move(read_step));
  }
  void add(assembly a) {
    if (m_rev_comp) {
      rev_asm(a);
    }
    m_cov->add(make_unique<assembly>(a));
  }
  read_id_set get_read_ids(const dna_sequence& seq) {
    seqset_range r = m_options.seqset->find(seq);
    CHECK(r.valid()) << seq;
    read_id_set read_ids;
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
  std::string print_read_ids(const read_id_set& read_ids) {
    std::stringstream os;

    if (read_ids.empty()) {
      return "(none)";
    }

    os << read_ids.size() << " reads: ";
    for (const auto& read_id : read_ids) {
      os << "\n" << read_id << "(";
      os << m_options.readmap->get_read_by_id(read_id).get_seqset_entry().sequence() << ")";
    }

    return os.str();
  }
  std::string print_all_read_ids(const assembly& a) {
    std::stringstream os;
    if (!a.edge_coverage) {
      return "(no edge coverage)";
    }
    const auto& ec = *a.edge_coverage;
    os << " var_start=" << print_read_ids(ec.variant_start);
    os << " var_end=" << print_read_ids(ec.variant_end);
    os << " ref_start=" << print_read_ids(ec.reference_start);
    os << " ref_end=" << print_read_ids(ec.reference_end);
    os << " interior=" << print_read_ids(ec.interior);
    return os.str();
  }
  void flush() {
    m_cov.reset();
    expect_sorted(assembly::left_offset_less_than);

    if (m_rev_comp) {
      reverse_found_assemblies();
      m_scaffold = m_scaffold.rev_comp();
    }

    for (const auto& a : m_ref_assemblies) {
      CHECK(a.edge_coverage) << a;
      m_ref_read_ids.insert(a.edge_coverage->interior.begin(), a.edge_coverage->interior.end());
    }
  }

  void rev_asm(assembly& a) {
    reverse_assembly_in_place(&a, m_options.readmap, m_options.scaffold->end_pos());
  }

 protected:
  bool m_rev_comp;
  scaffold m_rev_scaffold;
  boost::optional<add_ref> m_cov;
  aoffset_t m_ref_end_pos;
  read_id_set m_ref_read_ids;
};

TEST_P(pair_edge_cov_test, var_start) {
  use_pair_edge_ref({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdEFGH"), tseq_rc("uvwxyz")}});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcd").size();
  a.seq = tseq("EFGHI");
  a.right_offset = tseq("abcdefghi").size();
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmEdgeCoverIs(1, 0, 0, 0)));
}

TEST_P(pair_edge_cov_test, var_end) {
  use_pair_edge_ref({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("BCDEfgh"), tseq_rc("uvwxyz")}});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("ABCDE");
  a.right_offset = tseq("ABCDE").size();
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmEdgeCoverIs(0, 1, 0, 0)));
}

TEST_P(pair_edge_cov_test, ref_start) {
  use_pair_edge_ref({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdefgh"), tseq_rc("uvwxyz")}});

  start_calc();
  assembly a;
  a.left_offset = tseq("abcd").size();
  a.seq = tseq("EFGHI");
  a.right_offset = tseq("abcdefghi").size();
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmEdgeCoverIs(0, 0, 1, 0)));
}

TEST_P(pair_edge_cov_test, ref_end) {
  use_pair_edge_ref({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdefgh"), tseq_rc("uvwxyz")}});

  start_calc();
  assembly a;
  a.left_offset = 0;
  a.seq = tseq("ABCDE");
  a.right_offset = tseq("abcde").size();
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmEdgeCoverIs(0, 0, 0, 1)));
}

TEST_P(pair_edge_cov_test, insert) {
  use_pair_edge_ref({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdefgh"), tseq_rc("uvwxyz")},
                    {tseq("abcde") + dna_T + tseq("fghi"), tseq_rc("tuvwxy")}});

  auto ref_read_id = get_read_id(tseq("bcdefgh"));
  auto var_read_id = get_read_id(tseq("abcde") + dna_T + tseq("fghi"));

  start_calc();
  assembly a;
  a.left_offset = tseq("abcde").size();
  a.seq = dna_T;
  a.right_offset = tseq("abcde").size();
  add(a);
  flush();

  ASSERT_THAT(m_non_ref_assemblies, SizeIs(1));
  const auto& a_out = m_non_ref_assemblies[0];
  const auto& ec = *a_out.edge_coverage;
  EXPECT_THAT(ec.variant_start, UnorderedElementsAre(var_read_id)) << print_all_read_ids(a_out);
  EXPECT_THAT(ec.variant_end, UnorderedElementsAre(var_read_id)) << print_all_read_ids(a_out);
  EXPECT_THAT(ec.reference_start, UnorderedElementsAre(ref_read_id)) << print_all_read_ids(a_out);
  EXPECT_THAT(ec.reference_end, UnorderedElementsAre(ref_read_id)) << print_all_read_ids(a_out);
  EXPECT_THAT(ec.interior, IsEmpty());
}

// Test that we get all the proper read ids, and that they're
// base-correct and don't have any off by 1 errors.
TEST_P(pair_edge_cov_test, proper_read_ids) {
  use_pair_edge_ref(
      {{0,  //
        tseq("a") + dna_A + dna_T + tseq("bcdefgh") + dna_T + dna_A + tseq("ijklmnopqrstuvwxyz")}});
  dna_sequence ref_interior_seq = dna_T + tseq("bcdefgh") + dna_T;
  dna_sequence var_interior_seq = dna_G + tseq("BCDEFGH") + dna_G;
  dna_sequence ref_left_seq = dna_A + ref_interior_seq;
  dna_sequence var_left_seq = dna_A + var_interior_seq;
  dna_sequence ref_right_seq = ref_interior_seq + dna_A;
  dna_sequence var_right_seq = var_interior_seq + dna_A;

  dna_sequence outside_left_seq = tseq("a") + dna_A;
  dna_sequence ref_left_seq2 = tseq("a") + dna_A + dna_T;
  dna_sequence var_left_seq2 = tseq("a") + dna_A + dna_G;
  dna_sequence outside_right_seq = dna_A + tseq("ijklmn");
  dna_sequence ref_right_seq2 = dna_T + dna_A + tseq("ijkl");
  dna_sequence var_right_seq2 = dna_G + dna_A + tseq("ijkl");

  dna_sequence pair_support_seq = tseq_rc("uvwxyz");

  std::cout << "ref_left_seq: " << ref_left_seq << "\n";
  std::cout << "ref_left_seq2: " << ref_left_seq2 << "\n";
  std::cout << "ref_right_seq: " << ref_right_seq << "\n";
  std::cout << "ref_right_seq2: " << ref_right_seq2 << "\n";
  std::cout << "ref_interior_seq: " << ref_interior_seq << "\n";
  std::cout << "var_left_seq: " << var_left_seq << "\n";
  std::cout << "var_left_seq2: " << var_left_seq2 << "\n";
  std::cout << "var_right_seq: " << var_right_seq << "\n";
  std::cout << "var_right_seq2: " << var_right_seq2 << "\n";
  std::cout << "var_interior_seq: " << var_interior_seq << "\n";
  std::cout << "outside_left_seq: " << outside_left_seq << "\n";
  std::cout << "outside_right_seq: " << outside_right_seq << "\n";
  std::cout << "pair_support_seq: " << outside_right_seq << "\n";

  use_paired_reads({
      {ref_left_seq, pair_support_seq},      //
      {ref_left_seq2, pair_support_seq},     //
      {ref_right_seq, pair_support_seq},     //
      {ref_right_seq2, pair_support_seq},    //
      {ref_interior_seq, pair_support_seq},  //
      {var_left_seq, pair_support_seq},      //
      {var_left_seq2, pair_support_seq},     //
      {var_right_seq, pair_support_seq},     //
      {var_right_seq2, pair_support_seq},    //
      {var_interior_seq, pair_support_seq},  //
      {outside_left_seq, pair_support_seq},  //
      {outside_right_seq, pair_support_seq}  //
  });
  start_calc();
  assembly a;
  a.left_offset = tseq("a").size() + dna_T.size();
  a.seq = var_interior_seq;
  a.right_offset = tseq("a").size() + dna_T.size() + ref_interior_seq.size();
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmEdgeCoverIs(2, 2, 2, 2)));
  ASSERT_EQ(m_non_ref_assemblies.size(), 1);
  const auto& a_out = m_non_ref_assemblies[0];
  ASSERT_TRUE(a_out.edge_coverage);
  const auto& ec = *a_out.edge_coverage;
  EXPECT_THAT(ec.variant_start,
              UnorderedElementsAre(get_read_id(var_left_seq), get_read_id(var_left_seq2)))
      << print_all_read_ids(a_out);
  EXPECT_THAT(ec.variant_end,
              UnorderedElementsAre(get_read_id(var_right_seq), get_read_id(var_right_seq2)))
      << print_all_read_ids(a_out);
  EXPECT_THAT(ec.reference_start,
              UnorderedElementsAre(get_read_id(ref_left_seq), get_read_id(ref_left_seq2)))
      << print_all_read_ids(a_out);
  EXPECT_THAT(ec.reference_end,
              UnorderedElementsAre(get_read_id(ref_right_seq), get_read_id(ref_right_seq2)))
      << print_all_read_ids(a_out);
  EXPECT_THAT(ec.interior, UnorderedElementsAre(get_read_id(var_interior_seq)))  //
      << print_all_read_ids(a_out);

  read_id_set expected_ref_read_ids = get_read_ids(pair_support_seq.rev_comp());
  expected_ref_read_ids.insert(get_read_id(outside_left_seq));
  expected_ref_read_ids.insert(get_read_id(outside_right_seq));
  expected_ref_read_ids.insert(get_read_id(ref_interior_seq));

  EXPECT_THAT(m_ref_read_ids, UnorderedElementsAreArray(expected_ref_read_ids.to_vector()))
      << "Expecting\n"
      << print_read_ids(expected_ref_read_ids) << " in " << print_read_ids(m_ref_read_ids);
}

INSTANTIATE_TEST_CASE_P(pair_edge_cov_fwd_tests, pair_edge_cov_test, ::testing::Values(false));
INSTANTIATE_TEST_CASE_P(pair_edge_cov_rev_tests, pair_edge_cov_test, ::testing::Values(true));

class pair_edge_wild_test : public assemble_test {
 public:
  void start_calc() {
    auto edge_step = make_unique<pair_edge_cov>(m_options, test_output());
    auto pair_step = make_unique<pair_cov>(m_options, std::move(edge_step));
    auto read_step = make_unique<read_cov>(m_options, std::move(pair_step));
    m_cov.emplace(m_options, m_options.max_pair_distance, false /* whole ref */,
                  0 /* length limit */, std::move(read_step));
  }
  void add(assembly a) { m_cov->add(make_unique<assembly>(a)); }
  void flush() {
    m_cov.reset();
    expect_sorted(assembly::left_offset_less_than);
  }

  void add_vcf_assembly(const std::string& vcf_offset /* one based */, const std::string& ref,
                        const std::string& alt, int assembly_id = 0) {
    assembly a;
    a.assembly_id = assembly_id;
    a.left_offset = std::stoi(vcf_offset) - 1;
    a.right_offset = a.left_offset + ref.size();
    a.seq = alt;
    ASSERT_EQ(ref, m_scaffold.subscaffold_str(a.left_offset, a.right_offset - a.left_offset));

    // Trim anchor base if present
    if (a.seq.size() > 0 && !ref.empty() && a.seq.subseq(0, 1).as_string() == ref.substr(0, 1)) {
      a.seq = a.seq.subseq(1, a.seq.size() - 1);
      ++a.left_offset;

      ASSERT_EQ(ref.substr(1),
                m_scaffold.subscaffold_str(a.left_offset, a.right_offset - a.left_offset));
    }

    add(a);
  }

 protected:
  boost::optional<add_ref> m_cov;
};

TEST_F(pair_edge_wild_test, lambda_deletion_coverage) {
  size_t k_asm_id = 1;

  use_biograph("datasets/lambdaToyData/benchmark/father_lambda.bg");
  use_reference("datasets/lambdaToyData/benchmark/ref_lambda", "lambda");

  m_options.min_pair_distance = 0;
  m_options.max_pair_distance = 10000;

  start_calc();
  add_vcf_assembly("2191",
                   "TCTACGGAAAGCCGGTGGCCAGCATGCCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATAC"
                   "CGCGAAAGAGCAGATTTATAAC",
                   "T", k_asm_id);
  add_vcf_assembly("2667", "C", "CA");
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              UnorderedElementsAre(AsmEdgeCoverAndIdIs(121, 121, 0, 0, k_asm_id), AssemblyIdIs(0)));
}

}  // namespace variants

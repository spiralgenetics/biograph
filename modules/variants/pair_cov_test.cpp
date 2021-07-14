#include "modules/variants/pair_cov.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/place_pair_cov.h"
#include "modules/variants/read_cov.h"
#include "modules/variants/sort.h"

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
  CHECK(a.pair_read_coverage);
  return a.pair_read_coverage->reads();
}

class pair_cov_test : public assemble_test,
                      public WithParamInterface<std::pair<bool,  // rev_comp
                                                          bool   // Use pair placer
                                                          >> {
 public:
  pair_cov_test() {
    m_rev_comp = GetParam().first;
    m_pair_placer = GetParam().second;
  }

  Matcher<const assembly&> AsmPairReadCoverageIs(
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

  void start_calc() {
    m_ref_end_pos = m_options.scaffold->end_pos();

    if (m_rev_comp) {
      m_scaffold = m_scaffold.rev_comp();
      m_options.scaffold = &m_scaffold;
    }
    pipeline_step_t pcov;
    if (m_pair_placer) {
      m_popts.ideal_pair_distance = (m_options.max_pair_distance + m_options.min_pair_distance) / 2;
      pcov = make_unique<place_pair_cov>(m_options, m_popts, test_output());
    } else {
      pcov = make_unique<pair_cov>(m_options, test_output());
    }
    auto cov = make_unique<read_cov>(m_options, std::move(pcov));
    auto add = make_unique<add_ref>(m_options,
                                    m_options.max_pair_distance + m_options.seqset->max_read_len(),
                                    false /* whole ref */, 0 /* length limit */, std::move(cov));
    m_sorter.emplace(canon_assembly_order(), std::move(add));
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
  void pad(std::ostream& os, int n) {
    while (n > 0) {
      os << " ";
      --n;
    }
  }
  std::string print_read_cov(dna_slice seq, const boost::optional<read_coverage_t>& maybe_cov) {
    if (!maybe_cov) {
      return "(no coverage)";
    }
    std::stringstream os;

    const read_coverage_t& cov = *maybe_cov;

    if (cov.reads().empty()) {
      return "(empty coverage)";
    }

    aoffset_t npad = 0 - cov.reads().begin()->offset;
    if (npad < 0) {
      npad = 0;
    }
    os << cov.reads().size() << " reads:\n";
    pad(os, npad);
    os << seq.as_string() << "\n";
    for (const auto& rd : cov.reads()) {
      for (uint32_t read_id : rd.read_ids) {
        pad(os, npad + rd.offset);
        os << m_options.readmap->get_read_by_id(read_id).get_seqset_entry().sequence() << " (@"
           << rd.offset << ")\n";
      }
    }

    return os.str();
  }
  template <typename C>
  std::string print_all(const C& c) {
    std::stringstream os;
    for (const assembly& a : c) {
      os << "Assembly: " << a << "\nRead coverage: " << print_read_cov(a.seq, a.read_coverage)
         << "\nPair coverage: " << print_read_cov(a.seq, a.pair_read_coverage) << "\n";
    }
    return os.str();
  }

  std::string print_all_asms() {
    std::stringstream os;
    os << "\nREF assemblies:\n"
       << print_all(m_ref_assemblies) << "\nNON-ref assemblies:\n"
       << print_all(m_non_ref_assemblies);
    return os.str();
  }

  uint32_t get_read_id(const dna_sequence& seq) {
    auto ids = get_read_ids(seq);
    CHECK_EQ(ids.size(), 1) << seq;
    return *ids.begin();
  }
  void add(assembly a) {
    if (m_rev_comp) {
      rev_asm(a);
    }
    m_sorter->add(make_unique<assembly>(a));
  }
  void flush() {
    m_sorter.reset();
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
  }

  void rev_asm(assembly& a) {
    reverse_assembly_in_place(&a, m_options.readmap, m_options.scaffold->end_pos());
  }

 protected:
  bool m_rev_comp;
  bool m_pair_placer;
  boost::optional<sorter> m_sorter;
  aoffset_t m_ref_end_pos;
  place_pair_options m_popts;
};

TEST_P(pair_cov_test, simple) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = tseq("a").size();
  m_options.max_pair_distance = tseq("a").size() * 26;

  start_calc();
  assembly a;
  a.left_offset = tseq("a").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AsmPairReadCoverageIs(dna_G, {{dna_G + tseq("cdef"), 0}})))
      << print_all_asms();
}

TEST_P(pair_cov_test, insert) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = tseq("a").size();
  m_options.max_pair_distance = tseq("a").size() * 26;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AsmPairReadCoverageIs(dna_G, {{dna_G + tseq("cdef"), 0}})))
      << print_all_asms();
}

TEST_P(pair_cov_test, exceeds_max_distance) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() - 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() - 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmPairReadCoverageIs(dna_G, {})))
      << print_all_asms();
}

TEST_P(pair_cov_test, exceeds_max_distance_but_delete) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() - 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() - 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);

  assembly b;
  b.left_offset = tseq("abcdefg").size();
  b.right_offset = tseq("abcdefg").size() + 1;
  b.seq = dna_sequence();
  add(b);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AsmPairReadCoverageIs(dna_G, {{dna_G + tseq("cdef"), 0}}), _))
      << print_all_asms();
}

TEST_P(pair_cov_test, exceeds_min_distance) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() + 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() + 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmPairReadCoverageIs(dna_G, {})))
      << print_all_asms();
}

TEST_P(pair_cov_test, exceeds_min_distance_but_insert) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() + 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() + 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);

  assembly b;
  b.left_offset = tseq("abcdefg").size();
  b.right_offset = tseq("abcdefg").size();
  b.seq = dna_G;
  add(b);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AsmPairReadCoverageIs(dna_G, {{dna_G + tseq("cdef"), 0}}), _))
      << print_all_asms();
}

TEST_P(pair_cov_test, distance_ok) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size();
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size();

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AsmPairReadCoverageIs(dna_G, {{dna_G + tseq("cdef"), 0}})))
      << print_all_asms();
}

TEST_P(pair_cov_test, distance_exceeds_max) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() - 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() - 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmPairReadCoverageIs(dna_G, {})))
      << print_all_asms();
}

TEST_P(pair_cov_test, distance_exceeds_min) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), tseq_rc("lmnopq")}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() + 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq")).size() + 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmPairReadCoverageIs(dna_G, {})))
      << print_all_asms();
}

TEST_P(pair_cov_test, exceeds_max_distance_in_asm) {
  use_ref_parts({{0, tseq("abcd")}});
  use_paired_reads({{dna_G + tseq("cdef"), (tseq("lmnopq") + dna_G).rev_comp()}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq") + dna_G).size() - 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq") + dna_G).size() - 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G + tseq("cdefghijklmnopq") + dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmPairReadCoverageIs(a.seq, {})))
      << print_all_asms();
}

TEST_P(pair_cov_test, exceeds_min_distance_in_asm) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), (tseq("lmnopq") + dna_G).rev_comp()}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq") + dna_G).size() + 1;
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq") + dna_G).size() + 1;

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G + tseq("cdefghijklmnopq") + dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies, ElementsAre(AsmPairReadCoverageIs(a.seq, {})))
      << print_all_asms();
}

TEST_P(pair_cov_test, distance_ok_in_asm) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{dna_G + tseq("cdef"), (tseq("lmnopq") + dna_G).rev_comp()}});

  m_options.min_pair_distance = (dna_G + tseq("cdefghijklmnopq") + dna_G).size();
  m_options.max_pair_distance = (dna_G + tseq("cdefghijklmnopq") + dna_G).size();

  start_calc();
  assembly a;
  a.left_offset = tseq("ab").size();
  a.right_offset = tseq("ab").size();
  a.seq = dna_G + tseq("cdefghijklmnopq") + dna_G;
  add(a);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AsmPairReadCoverageIs(
                  a.seq, {{dna_G + tseq("cdef"), 0},
                          {(tseq("lmnopq") + dna_G), (dna_G + tseq("cdefghijk")).size()}})))
      << print_all_asms();
}

INSTANTIATE_TEST_CASE_P(pair_cov_fwd_tests, pair_cov_test,
                        ::testing::Values(std::make_pair(false, false)));
INSTANTIATE_TEST_CASE_P(pair_cov_rev_tests, pair_cov_test,
                        ::testing::Values(std::make_pair(true, false)));

INSTANTIATE_TEST_CASE_P(pair_placer_fwd_tests, pair_cov_test,
                        ::testing::Values(std::make_pair(false, true)));
INSTANTIATE_TEST_CASE_P(pair_placer_rev_tests, pair_cov_test,
                        ::testing::Values(std::make_pair(true, true)));

}  // namespace variants

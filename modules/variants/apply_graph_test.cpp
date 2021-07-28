#include "modules/variants/apply_graph.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/read_cov.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

// Copies of all assemblies so we can save them for test results
struct test_ctx {
  assembly a;
  read_coverage_t ref_coverage;
  std::vector<assembly> refs;

  optional_aoffset ref_left_offset;
  optional_aoffset ref_right_offset;

  scaffold ref_scaffold;
  edge_coverage_t edge_coverage;

  bool operator<(const test_ctx& rhs) const { return a.assembly_id < rhs.a.assembly_id; }
};

std::ostream& operator<<(std::ostream& os, const test_ctx& ctx) {
  return os << "\nContext assembly: " << PrintToString(ctx.a)
            << "\nCoverage: " << PrintToString(ctx.ref_coverage) << "\nRefs:\n"
            << PrintToString(ctx.refs);
}

Matcher<const std::vector<read_coverage_read_t>&> ReadCoverageIs(
    std::vector<read_coverage_read_t> cov) {
  return ContainerEq(cov);
}

std::set<std::pair<aoffset_t, aoffset_t>> deref_ref_bounds(const std::vector<assembly>& refs) {
  std::set<std::pair<aoffset_t, aoffset_t>> result;
  for (const auto& ref : refs) {
    CHECK(result.emplace(ref.left_offset, ref.right_offset).second)
        << "Duplicate reference bounds?";
  }
  return result;
}

aoffset_t deref_read_cov_len(const read_coverage_t& read_cov) { return read_cov.assembly_len(); }

const std::vector<read_coverage_read_t>& deref_read_cov(const read_coverage_t& read_cov) {
  return read_cov.reads();
}

const std::vector<read_coverage_read_t>& deref_optional_reads(
    const boost::optional<read_coverage_t>& read_cov) {
  return read_cov->reads();
}

class apply_graph_test : public assemble_test {
 public:
  Matcher<const test_ctx&> ContextIs(size_t aid, aoffset_t left_offset, aoffset_t right_offset,
                                     optional_aoffset ref_left_offset,
                                     optional_aoffset ref_right_offset,
                                     std::set<std::pair<dna_sequence, aoffset_t>> expected_covs,
                                     std::set<std::pair<aoffset_t, aoffset_t>> expected_bounds) {
    read_coverage_set expected;
    for (const auto& expected_seq : expected_covs) {
      const auto& seq = expected_seq.first;
      aoffset_t offset = expected_seq.second;

      for (uint32_t read_id : get_read_ids(seq)) {
        expected.insert(offset, read_id, m_options.readmap->get_readlength(read_id));
      }
    }

    read_coverage_t expected_cov = expected.build_and_clear(right_offset - left_offset);
    return AllOf(Field(&test_ctx::a, AllOf(Field(&assembly::assembly_id, Eq(aid)),
                                           Field(&assembly::left_offset, Eq(left_offset)),
                                           Field(&assembly::right_offset, Eq(right_offset)))),
                 Field(&test_ctx::ref_left_offset, Eq(ref_left_offset)),
                 Field(&test_ctx::ref_right_offset, Eq(ref_right_offset)),
                 Field(&test_ctx::ref_coverage,
                       AllOf(ResultOf(deref_read_cov_len, Eq(right_offset - left_offset)),
                             ResultOf(deref_read_cov, ContainerEq(expected_cov.reads())))),
                 Field(&test_ctx::refs, ResultOf(deref_ref_bounds, ContainerEq(expected_bounds))));
  }

  std::multiset<dna_sequence> deref_reads(const read_id_set& ids) {
    std::multiset<dna_sequence> result;

    for (uint32_t read_id : ids) {
      result.insert(m_options.readmap->get_read_by_id(read_id).get_seqset_entry().sequence());
    }

    return result;
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

  void start() {
    m_ctx = make_unique<apply_graph>(
        [this](graph_context ctx) {
          test_ctx tctx;

          tctx.a = *ctx.a;
          tctx.ref_coverage = ctx.ref_coverage();
          if (ctx.left_ref) {
            tctx.ref_left_offset = ctx.left_ref->left_offset;
          }
          if (ctx.right_ref) {
            tctx.ref_right_offset = ctx.right_ref->right_offset;
          }
          for (const auto& ref : ctx.refs) {
            tctx.refs.push_back(*ref);
          }
          tctx.ref_scaffold = ctx.ref_scaffold();
          tctx.ref_scaffold.save_all_storage();

          tctx.edge_coverage = ctx.differential_edge_coverage(
              tctx.ref_scaffold, *ctx.a->read_coverage, tctx.ref_coverage);

          m_actual.emplace_back(std::move(tctx));
        },
        test_output());
    m_read_cov = make_unique<read_cov>(m_options, std::move(m_ctx));
    m_add_ref = make_unique<add_ref>(m_options, 0, true /* whole ref */, 0, std::move(m_read_cov));
  }

  void add(size_t aid, aoffset_t left_offset, dna_sequence seq, aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = aid;
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    m_expected_non_ref.push_back(AssemblyIs(a->left_offset, a->seq, a->right_offset));
    m_add_ref->add(std::move(a));
  }

  void flush() {
    m_add_ref.reset();
    expect_sorted(assembly::left_offset_less_than);

    // None of this stuff should add or delete any assemblies.
    EXPECT_THAT(m_non_ref_assemblies, ElementsAreArray(m_expected_non_ref));

    std::sort(m_actual.begin(), m_actual.end());
  }

 protected:
  std::unique_ptr<apply_graph> m_ctx;
  std::unique_ptr<add_ref> m_add_ref;
  std::unique_ptr<read_cov> m_read_cov;

  std::vector<Matcher<const assembly&>> m_expected_non_ref;

  std::vector<test_ctx> m_actual;
};

// Copy extents into an array of pairs for easier matching
std::vector<std::pair<aoffset_t, dna_sequence>> get_extents(const scaffold& s) {
  std::vector<std::pair<aoffset_t, dna_sequence>> out;
  for (const auto& e : s.extents()) {
    out.push_back(std::make_pair(e.offset, dna_sequence(e.sequence.begin(), e.sequence.end())));
  }
  return out;
}

TEST_F(apply_graph_test, simple) {
  use_ref_parts({{0, tseq("abcdefghijklmnop")}});
  use_reads({tseq("abcdef"), tseq("efghi"), tseq("bc") + dna_T, dna_T + tseq("hi")});

  start();
  add(1, tseq("abc").size(), dna_T, tseq("abcdefg").size());
  flush();

  EXPECT_THAT(m_actual,
              ElementsAre(ContextIs(  // Aid, offsets:
                  1, tseq("abc").size(), tseq("abcdefg").size(),
                  // Surrounding ref offsets:
                  0, tseq("abcdefghijklmnop").size(),
                  // Ref coverage:
                  {{tseq("abcdef"), -tseq("abc").size()}, {tseq("efghi"), tseq("d").size()}},
                  // Ref bounds:
                  {{tseq("abc").size(), tseq("abcdefg").size()}})));

  auto& ctx = m_actual[0];
  scaffold s = ctx.ref_scaffold;
  EXPECT_EQ(s.end_pos(), tseq("defg").size());
  EXPECT_THAT(get_extents(s), ElementsAre(Pair(0, tseq("defg"))));

  edge_coverage_t ec = ctx.edge_coverage;
  EXPECT_EQ(0, ec.start_common);
  EXPECT_EQ(0, ec.end_common);

  EXPECT_THAT(deref_reads(ec.variant_start), ElementsAre(tseq("bc") + dna_T));
  EXPECT_THAT(deref_reads(ec.variant_end), ElementsAre(dna_T + tseq("hi")));

  EXPECT_THAT(deref_reads(ec.interior), IsEmpty());
  EXPECT_THAT(deref_reads(ec.reference_start), ElementsAre(tseq("abcdef")));
  EXPECT_THAT(deref_reads(ec.reference_end), ElementsAre(tseq("efghi")));
}

TEST_F(apply_graph_test, multi_with_decoy) {
  use_ref_parts({{0, tseq("abcdefgh")}, {100, tseq("ABCDEFGH")}});

  use_reads({
      // Should be included in ref coverage:
      tseq("def"),
      tseq("efg"),
      tseq("ABC"),
      tseq("BCD"),
      // Should not be included in ref coverage:
      tseq("bc") + dna_T,
      tseq("cd") + dna_T,
      dna_T + tseq("CD"),
      dna_T + tseq("DE"),
      tseq("c") + dna_T + tseq("D"),
      tseq("cd") + dna_T + tseq("CD"),
  });

  start();
  add(1, tseq("abc").size(), dna_T, 100 + tseq("ABC").size());
  add(2, tseq("abcd").size(), dna_T, 100 + tseq("AB").size());
  flush();

  aoffset_t gap_size = 100 - tseq("abcdefgh").size();

  // Make sure the non-ref reads actually got included in coverage.
  EXPECT_THAT(m_non_ref_assemblies, Each(Field(&assembly::read_coverage,
                                               ResultOf(deref_optional_reads, Not(IsEmpty())))));

  EXPECT_THAT(m_actual,
              ElementsAre(ContextIs(  // Aid, offsets:
                              1, tseq("abc").size(), 100 + tseq("ABC").size(),
                              // Surrounding ref offsets:
                              0, 100 + tseq("ABCDEFGH").size(),
                              // Ref coverage:
                              {{tseq("def"), 0},
                               {tseq("efg"), tseq("d").size()},
                               {tseq("ABC"), tseq("defgh").size() + gap_size},
                               {tseq("BCD"), tseq("defgh").size() + gap_size + tseq("A").size()}},
                              // Ref bounds:
                              {{tseq("abc").size(), tseq("abcd").size()},
                               {tseq("abcd").size(), tseq("abcdefgh").size()},
                               {100, 100 + tseq("AB").size()},
                               {100 + tseq("AB").size(), 100 + tseq("ABC").size()}}),
                          ContextIs(  // Aid, offsets:
                              2, tseq("abcd").size(), 100 + tseq("AB").size(),
                              // Surrounding ref offsets:
                              tseq("abc").size(), 100 + tseq("ABC").size(),
                              // Ref coverage:
                              {{tseq("def"), -tseq("d").size()},
                               {tseq("efg"), 0},
                               {tseq("ABC"), tseq("efgh").size() + gap_size},
                               {tseq("BCD"), tseq("efgh").size() + gap_size + tseq("A").size()}},
                              // Ref bounds:
                              {{tseq("abcd").size(), tseq("abcdefgh").size()},
                               {100, 100 + tseq("AB").size()}})));
}

TEST_F(apply_graph_test, insert) {
  use_ref_parts({{0, tseq("abcdefgh")}});

  use_reads({tseq("cd") + dna_T + tseq("ef"), tseq("defg")});

  start();
  add(1, tseq("abcd").size(), dna_T, tseq("abcd").size());
  flush();

  EXPECT_THAT(m_actual, ElementsAre(ContextIs(
                            // Aid, offsets:
                            1, tseq("abcd").size(), tseq("abcd").size(),
                            // Surrounding ref offsets:
                            0, tseq("abcdefgh").size(),
                            // Ref coverage:
                            {{tseq("defg"), -tseq("d").size()}},
                            // Ref bounds:
                            {})));
}

TEST_F(apply_graph_test, deletion) {
  use_ref_parts({{0, tseq("abcdefgh")}});

  use_reads({tseq("cdgh"), tseq("defg")});

  start();
  add(1, tseq("abcd").size(), dna_sequence(), tseq("abcdef").size());
  flush();

  EXPECT_THAT(m_actual, ElementsAre(ContextIs(
                            // Aid, offsets:
                            1, tseq("abcd").size(), tseq("abcdef").size(),
                            // Surrounding ref offsets:
                            0, tseq("abcdefgh").size(),
                            // Ref coverage:
                            {{tseq("defg"), -tseq("d").size()}},
                            // Ref bounds:
                            {{tseq("abcd").size(), tseq("abcdef").size()}})));
}

TEST_F(apply_graph_test, big_insert) {
  use_ref_parts({{0, tseq("abcdefgh")}});

  use_reads({
      // Var:
      tseq("bcd"),
      tseq("cd") + dna_T,
      dna_T + tseq("AB"),
      tseq("ABC"),
      tseq("CDE"),
      tseq("DE") + dna_T,
      dna_T + tseq("ef"),
      tseq("efg"),
      // Ref:
      tseq("bcdef"),
  });

  start();
  add(1, tseq("abcd").size(), dna_T + tseq("ABCDE") + dna_T, tseq("abcd").size());
  flush();

  EXPECT_THAT(m_actual, ElementsAre(ContextIs(
                            // Aid, offsets:
                            1, tseq("abcd").size(), tseq("abcd").size(),
                            // Surrounding ref offsets:
                            0, tseq("abcdefgh").size(),
                            // Ref coverage:
                            {{tseq("bcdef"), -tseq("bcd").size()}},
                            // Ref bounds:
                            {})));

  auto& ctx = m_actual[0];
  scaffold s = ctx.ref_scaffold;
  EXPECT_EQ(s.end_pos(), 0);

  edge_coverage_t ec = ctx.edge_coverage;
  EXPECT_EQ(0, ec.start_common);
  EXPECT_EQ(0, ec.end_common);

  EXPECT_THAT(deref_reads(ec.variant_start), ElementsAre(tseq("cd") + dna_T));
  EXPECT_THAT(deref_reads(ec.variant_end), ElementsAre(dna_T + tseq("ef")));

  EXPECT_THAT(deref_reads(ec.interior), UnorderedElementsAre(dna_T + tseq("AB"), tseq("ABC"),
                                                             tseq("CDE"), tseq("DE") + dna_T));

  EXPECT_THAT(deref_reads(ec.reference_start), ElementsAre(tseq("bcdef")));
  EXPECT_THAT(deref_reads(ec.reference_end), ElementsAre(tseq("bcdef")));
}

// TODO(nils): Test when ec.start_common and ec.end_common are nonzero

}  // namespace variants

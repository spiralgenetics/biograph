#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/variants/discovery/discovery_testutil.h"
#include "modules/variants/discovery/state.h"

#include <fstream>

using namespace testing;
using namespace dna_testutil;

namespace variants {
namespace discovery {

class walk_ref_test : public discovery_test, public WithParamInterface<bool> {
 public:
  walk_ref_test() { m_rev_comp = GetParam(); }
  void add_expected_ref(aoffset_t start, dna_sequence seq) {
    m_expected_ref.emplace(m_rev_comp, start, seq);
    m_expected_ref.emplace(!m_rev_comp,
                           m_options.scaffold->end_pos() - (start + aoffset_t(seq.size())),
                           seq.rev_comp());
  }

  void run_walk_ref_test() {
    init_discovery();
    add_ref();
    save_ref_locations();
    save_search_entries();
    save_pair_support();
    EXPECT_THAT(m_pop_entries, IsEmpty());
    EXPECT_THAT(m_rejoin_entries, IsEmpty());
  }
  void save_ref_locations() {
    for (view_t* v : m_st->both_dirs()) {
      for (const auto& r_and_ri : range_info_table(v)) {
        const seqset_range& r = r_and_ri.first;
        const range_info_t& ri = r_and_ri.second;

        for (aoffset_t ref_loc : ri.reference_offsets) {
          m_actual_ref.emplace(v->is_rev_comp(), ref_loc, r.sequence());
        }
      }
    }
  }

  std::set<std::tuple<bool /* is reverse */, aoffset_t, dna_sequence>> m_actual_ref;
  std::set<std::tuple<bool /* is reverse */, aoffset_t, dna_sequence>> m_expected_ref;
};

TEST_P(walk_ref_test, simple) {
  m_options.min_overlap = tseq("efgh").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("efghi"), tseq("ghijk")});
  add_expected_ref(tseq("abcd").size(), tseq("efghi"));
  add_expected_ref(tseq("abcdef").size(), tseq("ghijk"));
  run_walk_ref_test();

  EXPECT_THAT(m_push_entries, IsEmpty());
  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));
}

TEST_P(walk_ref_test, branch_fwd) {
  m_options.min_overlap = tseq("efgh").size();
  use_ref_parts({{0, tseq("abcdefg") + dna_A + tseq("hijklmnopqrstuvwxyz")}});
  use_reads({tseq("DEFG") + dna_T + tseq("hijk")});
  add_expected_ref(tseq("abcdefg").size() + 1, tseq("hijk"));
  run_walk_ref_test();

  EXPECT_THAT(m_push_entries, ElementsAre(FwdPushSearchEntry(
                                  tseq("hijk").size(), tseq("abcdefg").size() + dna_A.size(),
                                  dna_T + tseq("hijk"), dna_T + tseq("hijk"))));
  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));
}

TEST_P(walk_ref_test, branch_rev) {
  m_options.min_overlap = tseq("efgh").size();
  use_ref_parts({{0, tseq("abcdefg") + dna_A + tseq("hijklmnopqrstuvwxyz")}});
  use_reads({tseq("defg") + dna_T + tseq("HIJK"), tseq("defg")});
  add_expected_ref(tseq("abc").size(), tseq("defg"));
  run_walk_ref_test();

  EXPECT_THAT(m_push_entries,
              ElementsAre(RevPushSearchEntry(tseq("defg").size(), tseq("abcdefg").size(),
                                             tseq("defg") + dna_T, tseq("defg") + dna_T)));
  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));
}

TEST_P(walk_ref_test, not_enough_overlap_edges) {
  // Make minimum overlap one base longer.
  m_options.min_overlap = tseq("efgh").size() + 1;
  use_ref_parts({{0, tseq("abcdefg") + dna_A + tseq("hijklmnopqrstuvwxyz")}});
  use_reads({tseq("defg") + dna_T + tseq("HIJK"), tseq("DEFG") + dna_T + tseq("hijk")});
  run_walk_ref_test();
  EXPECT_THAT(m_push_entries, IsEmpty());
  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));
}

TEST_P(walk_ref_test, multi_extent_edges_exactly_min_overlap) {
  m_options.min_overlap = tseq("abcd").size();

  use_ref_parts(
      {{1000, tseq("abcdefghijklmnopqrstuvwxyz")}, {2000, tseq("ABCDEFGHIJKLMNOPQRSTUVWXYZ")}});

  add_expected_ref(1000, tseq("abcd"));
  add_expected_ref(1000 + tseq("abcdefghijklmnopqrstuv").size(), tseq("wxyz"));
  add_expected_ref(2000, tseq("ABCD"));
  add_expected_ref(2000 + tseq("ABCDEFGHIJKLMNOPQRSTUV").size(), tseq("WXYZ"));

  use_reads({tseq("abcd"), tseq("ABCD"), tseq("wxyz"), tseq("WXYZ")});
  run_walk_ref_test();

  EXPECT_THAT(m_push_entries, IsEmpty());
  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));
}

TEST_P(walk_ref_test, multi_extent_edges_extend) {
  m_options.min_overlap = tseq("abcd").size();

  use_ref_parts(
      {{1000, tseq("abcdefghijklmnopqrstuvwxyz")}, {2000, tseq("ABCDEFGHIJKLMNOPQRSTUVWXYZ")}});

  add_expected_ref(1000, tseq("abcd"));
  add_expected_ref(1000 + tseq("abcdefghijklmnopqrstuv").size(), tseq("wxyz"));
  add_expected_ref(2000, tseq("ABCD"));
  add_expected_ref(2000 + tseq("ABCDEFGHIJKLMNOPQRSTUV").size(), tseq("WXYZ"));

  use_reads(
      {dna_A + tseq("abcd"), dna_A + tseq("ABCD"), tseq("wxyz") + dna_A, tseq("WXYZ") + dna_A});
  run_walk_ref_test();

  EXPECT_THAT(
      m_push_entries,
      ElementsAre(
          FwdPushSearchEntry(tseq("abcd").size(), 1000, dna_A + tseq("abcd"), dna_A + tseq("abcd")),
          FwdPushSearchEntry(tseq("ABCD").size(), 2000, dna_A + tseq("ABCD"), dna_A + tseq("ABCD")),
          RevPushSearchEntry(tseq("WXYZ").size(), 2000 + tseq("ABCDEFGHIJKLMNOPQRSTUVWXYZ").size(),
                             tseq("WXYZ") + dna_A, tseq("WXYZ") + dna_A),
          RevPushSearchEntry(tseq("wxyz").size(), 1000 + tseq("abcdefghijklmnopqrstuvwxyz").size(),
                             tseq("wxyz") + dna_A, tseq("wxyz") + dna_A)));

  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));
}

TEST_P(walk_ref_test, multi_extent_edges_less_than_min_overlap) {
  m_options.min_overlap = tseq("abcd").size() + 1;

  use_ref_parts(
      {{1000, tseq("abcdefghijklmnopqrstuvwxyz")}, {2000, tseq("ABCDEFGHIJKLMNOPQRSTUVWXYZ")}});

  use_reads({tseq("abcd"), tseq("ABCD"), tseq("wxyz"), tseq("WXYZ")});

  run_walk_ref_test();
  EXPECT_THAT(m_push_entries, IsEmpty());
  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));
}

TEST_P(walk_ref_test, pair_position) {
  m_options.min_overlap = tseq("abcd").size();
  m_options.min_pair_distance = 100;
  m_options.max_pair_distance = 200;

  use_ref_parts({{1000, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("abcde"), tseq_rc("1234")}, {tseq("56789"), tseq_rc("wxyz")}}, {});

  add_expected_ref(1000, tseq("abcde"));
  add_expected_ref(1000 + tseq("abcdefghijklmnopqrstuv").size(), tseq("wxyz"));
  run_walk_ref_test();

  EXPECT_THAT(m_push_entries, IsEmpty());
  //  EXPECT_THAT(m_actual_ref, ContainerEq(m_expected_ref));

  EXPECT_THAT(m_pair_support, SizeIs(2));
  interval_set_t expected_1234;
  expected_1234 +=
      interval_t(1000 + 100 - tseq_rc("1234").size(), 1000 + 200 - tseq("1234").size());
  EXPECT_EQ(m_pair_support[tseq("1234")], expected_1234);

  interval_set_t expected_56789;
  expected_56789 += interval_t(1000 + tseq("abcdefghijklmnopwrstuvwxyz").size() - 200,
                               1000 + tseq("abcdefghijklmnopwrstuvwxyz").size() - 100);
  EXPECT_EQ(m_pair_support[tseq("56789")], expected_56789);
}

INSTANTIATE_TEST_CASE_P(fwd_walk_ref_test, walk_ref_test,
                        ::testing::Values(false /* not rev_comp */));
INSTANTIATE_TEST_CASE_P(rev_walk_ref_test, walk_ref_test, ::testing::Values(true /* rev_comp */));

}  // namespace discovery
}  // namespace variants

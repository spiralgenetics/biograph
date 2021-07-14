#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/variants/discovery/discovery_testutil.h"
#include "modules/variants/discovery/state.h"
#include "modules/variants/discovery/state.h"

#include <fstream>

using namespace testing;
using namespace dna_testutil;

namespace variants {
namespace discovery {

class pop_search_test : public discovery_test, public WithParamInterface<bool /* rev comp */> {
 public:
  pop_search_test() { m_rev_comp = GetParam(); }

  void run_fwd_search_entry(int anchor_len, aoffset_t left_offset, dna_slice seq, dna_slice r) {
    init_discovery();
    add_ref_without_search();
    {
      CHECK_GT(seq.size(), anchor_len);
      std::cout << "Seq " << seq << "\n";
      seqset_range rc_r = m_options.seqset->find(seq.subseq(0, anchor_len).rev_comp());
      CHECK(rc_r.valid()) << seq.subseq(0, anchor_len);

      path rc_path(m_st->opts().readmap, rc_r.sequence(), rc_r, anchor_len,
                   0 /* bases since read */, anchor_len);

      std::cout << "RC path so far: " << rc_path << "\n";

      rc_path.push_front_drop(seq.subseq(anchor_len, seq.size() - anchor_len).rev_comp());
      std::cout << "Now rc path is: " << rc_path << "\n";

      branch* br = rev_view()->get_branch(seq[anchor_len].complement(),
                                          rev_view()->reverse_offset(left_offset));
      br->check_path_invariants(rc_path);

      seqset_range popped = get_seqset_range(r);
      CHECK_EQ(popped.sequence(), rc_path.range().sequence().subseq(0, popped.size()).rev_comp());

      std::unique_ptr<pop_search_entry> e =
          make_unique<pop_search_entry>(popped, rc_path, 0 /* no pair matches */);
      e->check_invariants(br);
      execute_search(br, std::move(e));
      br->check_invariants();
    }
    m_st->check_invariants();
    save_search_entries();
    save_partials();
  }
};

TEST_P(pop_search_test, simple_rejoin_ref) {
  m_options.min_overlap = m_options.min_pop_overlap = tseq("efgh").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("bcde") + dna_T + tseq("efgh"), tseq("efghi")});

  run_fwd_search_entry(tseq("bcde").size(), tseq("abcde").size(),
                       tseq("bcde") + dna_T + tseq("efghi"), tseq("efghi"));

  EXPECT_THAT(m_right_partials, IsEmpty());
  EXPECT_THAT(m_left_partials, IsEmpty());
  EXPECT_THAT(m_push_entries, IsEmpty());
  // EXPECT_THAT(m_pop_entries, IsEmpty());
  EXPECT_THAT(m_rejoin_entries, ElementsAre(RejoinSearchEntry(
                                    tseq("efgh").size(), tseq("abcde").size(),
                                    tseq("bcde") + dna_T + tseq("efghi"), tseq("abcd").size())));
}

TEST_P(pop_search_test, rejoin_ref_after_pop) {
  m_options.min_overlap = m_options.min_pop_overlap = tseq("efg").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("bcde") + dna_T + dna_A + tseq("efgh"), tseq("efghi")});

  dna_sequence orig_seq = tseq("bcde") + dna_T + dna_A + tseq("efg");
  unsigned path_overlap = tseq("bcde").size();
  dna_sequence r_seq = dna_T + dna_A + tseq("efg");
  run_fwd_search_entry(path_overlap, tseq("abcde").size(), orig_seq, r_seq);

  EXPECT_THAT(m_pop_entries,
              ElementsAre(FwdPopSearchEntry(r_seq.size() - 1, tseq("abcde").size(), orig_seq,
                                            r_seq.subseq(1, r_seq.size() - 1))));

  search_each_branch_once();
  save_search_entries();

  EXPECT_THAT(m_pop_entries,
              ElementsAre(FwdPopSearchEntry(r_seq.size() - 2, tseq("abcde").size(), orig_seq,
                                            r_seq.subseq(2, r_seq.size() - 2))));

  search_each_branch_once();
  save_search_entries();

  EXPECT_THAT(m_right_partials, IsEmpty());
  EXPECT_THAT(m_left_partials, IsEmpty());
  EXPECT_THAT(m_push_entries, IsEmpty());
  //  EXPECT_THAT(m_pop_entries, IsEmpty());
  EXPECT_THAT(m_rejoin_entries,
              ElementsAre(RejoinSearchEntry(tseq("efg").size(), tseq("abcde").size(),
                                            tseq("bcde") + dna_T + dna_A + tseq("efg"),
                                            tseq("abcd").size())));
}

TEST_P(pop_search_test, convert_to_push) {
  // A pop entry that finds pair support should convert to a push search.

  m_options.min_overlap = tseq("1234").size();
  m_options.min_pop_overlap = tseq("12").size();
  m_options.min_pair_distance = 100;
  m_options.max_pair_distance = 2000;

  use_ref_parts(
      {{0, tseq("abcdefghijklmnopqrstuvwxyz")}, {1000, tseq("ABCDEFGHIJKLMNOPQRSTUVWXYZ")}});
  dna_sequence orig_seq = tseq("bcdef") + dna_T + dna_A + tseq("12");
  dna_sequence r_seq = dna_T + dna_A + tseq("12");
  use_paired_reads({{tseq("12345"), tseq_rc("ABCD")}},
                   {orig_seq,
                    // Intermediate read that should cause the total path overlap not to get smaller
                    // than tseq("123").size()
                    tseq("def") + dna_T + dna_A + tseq("123")});

  run_fwd_search_entry(tseq("bcdef").size(), tseq("abcdef").size(), orig_seq, r_seq);

  EXPECT_THAT(m_pop_entries,
              ElementsAre(FwdPopSearchEntry(r_seq.size() - 1, tseq("abcdef").size(), orig_seq,
                                            r_seq.subseq(1, r_seq.size() - 1))));

  search_each_branch_once();
  save_search_entries();

  EXPECT_THAT(m_pop_entries,
              ElementsAre(FwdPopSearchEntry(r_seq.size() - 2, tseq("abcdef").size(), orig_seq,
                                            r_seq.subseq(2, r_seq.size() - 2))));

  search_each_branch_once();
  save_search_entries();

  EXPECT_THAT(m_right_partials, IsEmpty());
  EXPECT_THAT(m_left_partials, IsEmpty());
  EXPECT_THAT(m_push_entries, ElementsAre(RevPushSearchEntry(
                                  tseq("123").size(), tseq("abcdef").size(),
                                  tseq("bcdef") + dna_T + dna_A + tseq("12345"), tseq("12345"))));
  EXPECT_THAT(m_pop_entries, IsEmpty());
  EXPECT_THAT(m_rejoin_entries, IsEmpty());
}

INSTANTIATE_TEST_CASE_P(fwd_pop_search_test, pop_search_test,
                        ::testing::Values(false /* not rev_comp */));
INSTANTIATE_TEST_CASE_P(rev_pop_search_test, pop_search_test,
                        ::testing::Values(true /* rev_comp */));

}  // namespace discovery
}  // namespace variants

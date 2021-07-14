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

class push_search_test : public discovery_test, public WithParamInterface<bool /* rev comp */> {
 public:
  push_search_test() {
    m_rev_comp = GetParam();
    // Make sure we don't use the pop tracer by accident.
    m_options.min_pop_overlap = 1000;
  }

  void init_push() {
    init_discovery();
    add_ref_without_search();
  }

  void run_fwd_search_entry(unsigned anchor_size, aoffset_t right_offset, dna_sequence seq,
                            dna_sequence r) {
    init_push();
    add_fwd_search_entry(anchor_size, right_offset, seq, r);
  }

  void add_fwd_search_entry(unsigned anchor_size, aoffset_t right_offset, dna_sequence seq,
                            dna_sequence r) {
    view_t* v = fwd_view();
    path p(m_st->opts().readmap, seq, get_seqset_range(r), anchor_size, 0 /* bases since read */,
           anchor_size);

    branch* br = v->get_branch(seq[seq.size() - anchor_size - 1], right_offset - anchor_size);

    std::unique_ptr<push_search_entry> e =
        make_unique<push_search_entry>(std::move(p), 0 /* no pair matches */);
    e->check_invariants(br);
    execute_search(br, std::move(e));
    br->check_invariants();
    m_st->check_invariants();
    save_search_entries();
    save_partials();
  }
  void add_rev_search_entry(unsigned anchor_size, aoffset_t left_offset, dna_sequence seq,
                            dna_sequence r) {
    view_t* v = rev_view();

    branch* br =
        v->get_branch(seq[anchor_size].complement(), v->reverse_offset(left_offset + anchor_size));

    path p(m_st->opts().readmap, seq.rev_comp(), get_seqset_range(r.rev_comp()), anchor_size,
           0 /* bases since read */, anchor_size);

    std::unique_ptr<push_search_entry> e =
        make_unique<push_search_entry>(std::move(p), 0 /* no pair matches */);
    e->check_invariants(br);
    execute_search(br, std::move(e));
    br->check_invariants();
    m_st->check_invariants();
    save_search_entries();
    save_partials();
  }
};

TEST_P(push_search_test, simple) {
  m_options.min_overlap = tseq("efgh").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcde"), tseq("bcde") + dna_T + tseq("efgh")});
  run_fwd_search_entry(tseq("efgh").size(), tseq("abcdefgh").size(), dna_T + tseq("efgh"),
                       dna_T + tseq("efgh"));

  //  EXPECT_THAT(m_right_partials, IsEmpty());
  EXPECT_THAT(m_left_partials, IsEmpty());
  // EXPECT_THAT(m_push_entries, IsEmpty());
  //  EXPECT_THAT(m_pop_entries, IsEmpty());
  EXPECT_THAT(m_rejoin_entries, ElementsAre(RejoinSearchEntry(
                                    tseq("efgh").size(), tseq("abcde").size(),
                                    tseq("abcde") + dna_T + tseq("efgh"), tseq("abcd").size())));
}

TEST_P(push_search_test, hanging_end) {
  m_options.min_overlap = tseq("efgh").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("bcde") + dna_T + tseq("efgh")});
  run_fwd_search_entry(tseq("efgh").size(), tseq("abcdefgh").size(), dna_T + tseq("efgh"),
                       dna_T + tseq("efgh"));

  EXPECT_THAT(m_right_partials,
              ElementsAre(Pair(tseq("abcdefgh").size(), tseq("bcde") + dna_T + tseq("efgh"))));
  EXPECT_THAT(m_left_partials, IsEmpty());
  EXPECT_THAT(m_push_entries, IsEmpty());
  EXPECT_THAT(m_pop_entries, ElementsAre(RevPopSearchEntry(tseq("efgh").size(), tseq("abcd").size(),
                                                           tseq("bcde") + dna_T + tseq("efgh"),
                                                           tseq("bcde") + dna_T + tseq("efgh"))));
  EXPECT_THAT(m_rejoin_entries, IsEmpty());
}

TEST_P(push_search_test, join_dead_ends_via_pop) {
  m_options.min_overlap = tseq("ghij").size();
  m_options.min_pop_overlap = tseq("12").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("bcde") + dna_T + tseq("12"), tseq("12") + dna_T + tseq("ghij")});
  init_push();
  add_fwd_search_entry(tseq("bcde").size(), tseq("abcdefghij").size(), dna_T + tseq("ghij"),
                       dna_T + tseq("ghij"));
  add_rev_search_entry(tseq("ghij").size(), tseq("a").size(), tseq("bcde") + dna_T,
                       tseq("bcde") + dna_T);

  save_search_entries();
  save_partials();
  EXPECT_THAT(m_pop_entries, SizeIs(2));
  while (m_rejoin_entries.empty() && !m_pop_entries.empty()) {
    search_each_branch_once();
    save_search_entries();
  }
  EXPECT_THAT(m_push_entries, IsEmpty());
  // EXPECT_THAT(m_pop_entries, IsEmpty());

  // Both branches should generate rejoins.
  auto rejoin_matcher = RejoinSearchEntry(tseq("12").size(), tseq("abcde").size(),
                                          tseq("bcde") + dna_T + tseq("12") + dna_T + tseq("ghij"),
                                          tseq("abcdef").size());
  EXPECT_THAT(m_rejoin_entries, ElementsAre(rejoin_matcher, rejoin_matcher));
}

TEST_P(push_search_test, join_right_partial) {
  m_options.min_overlap = tseq("gh").size();
  m_options.ref_align_factor = 10;
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("e") + dna_T + tseq("12"), tseq("1234"), tseq("34") + dna_A + tseq("ghij")});

  init_push();
  view_t* v = fwd_view();
  right_partial rc_rp((tseq("e") + dna_T + tseq("1234")).rev_comp(),
                      v->reverse_offset(tseq("abcd").size()), 0 /* no pair matches */);

  v->reverse_view()->add_right_partial(get_seqset_range(tseq_rc("1234")), rc_rp);

  add_fwd_search_entry(tseq("ghij").size(), tseq("abcdefghij").size(), dna_A + tseq("ghij"),
                       dna_A + tseq("ghij"));

  // Should have advanced before it had to decrease overlap.
  EXPECT_THAT(m_push_entries,
              ElementsAre(FwdPushSearchEntry(
                  tseq("34").size(), tseq("abcdef").size(),
                  drop_front(tseq("12").size() - 1, tseq("1234") + dna_A + tseq("ghij")),
                  drop_front(tseq("12").size() - 1, tseq("1234")))));
  // EXPECT_THAT(m_pop_entries, IsEmpty());
  EXPECT_THAT(m_rejoin_entries, IsEmpty());

  // Continue trying to find the right partial in the other direction
  while (m_rejoin_entries.empty() && !m_push_entries.empty()) {
    search_each_branch_once();
    save_search_entries();
  }
  // EXPECT_THAT(m_push_entries, IsEmpty());
  // EXPECT_THAT(m_pop_entries, IsEmpty());
  EXPECT_THAT(m_rejoin_entries,
              ElementsAre(RejoinSearchEntry(tseq("e").size(), tseq("abcde").size(),
                                            tseq("e") + dna_T + tseq("1234") + dna_A + tseq("ghij"),
                                            tseq("abcdef").size())));
}

INSTANTIATE_TEST_CASE_P(fwd_push_search_test, push_search_test,
                        ::testing::Values(false /* not rev_comp */));
INSTANTIATE_TEST_CASE_P(rev_push_search_test, push_search_test,
                        ::testing::Values(true /* rev_comp */));

}  // namespace discovery
}  // namespace variants

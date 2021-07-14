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

class rejoin_test : public discovery_test, public WithParamInterface<bool /* rev comp */> {
 public:
  rejoin_test() { m_rev_comp = GetParam(); }

  void run_fwd_search_entry(int path_overlap, aoffset_t left_offset, aoffset_t left_anchor_len,
                            aoffset_t right_offset, aoffset_t right_anchor_len, dna_slice seq) {
    init_discovery();
    add_ref_without_search();
    {
      dna_slice right_anchor = seq.subseq(seq.size() - right_anchor_len, right_anchor_len);
      seqset_range init_r = m_options.seqset->find(right_anchor);
      CHECK(init_r.valid()) << right_anchor;
      path p(m_options.readmap, right_anchor, init_r, right_anchor_len, 0, right_anchor_len);
      p.push_front_drop(seq.subseq(0, seq.size() - right_anchor_len));
      branch* br = fwd_view()->get_branch(seq[seq.size() - right_anchor_len - 1], right_offset);
      std::unique_ptr<rejoin_search_entry> e = make_unique<rejoin_search_entry>(
          path_overlap, left_offset, left_anchor_len, std::move(p), 0 /* no pair matches */);
      e->check_invariants(br);
      execute_search(br, std::move(e));
      br->check_invariants();
    }
    m_st->check_invariants();
    save_search_entries();
    save_partials();
    EXPECT_THAT(m_right_partials, IsEmpty());
    EXPECT_THAT(m_left_partials, IsEmpty());
    EXPECT_THAT(m_pop_entries, IsEmpty());
    EXPECT_THAT(m_rejoin_entries, IsEmpty());

    if (m_rev_comp) {
      reverse_found_assemblies();
    }
  }
};

TEST_P(rejoin_test, simple_rejoin_ref) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("bcde") + dna_T + tseq("ghi")});

  run_fwd_search_entry(tseq("ghi").size(), tseq("abcde").size(), tseq("cde").size(),
                       tseq("abcdef").size(), tseq("ghi").size(),
                       tseq("cde") + dna_T + tseq("ghi"));

  EXPECT_THAT(m_push_entries, IsEmpty());

  EXPECT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(tseq("ab").size(), tseq("cde") + dna_T + tseq("ghi"),
                                     tseq("abcdefghi").size())));
}

TEST_P(rejoin_test, rejoin_and_search_more) {
  m_options.min_overlap = tseq("ghi").size();
  use_ref_parts({{0, tseq("abcde") + dna_T + tseq("fghijklmnopqrstuvwxyz")}});
  aoffset_t f_pos = (tseq("abcde") + dna_T).size();
  use_reads({
      tseq("bcde") + dna_A + tseq("ghij"),  // rejoining variant
      dna_G + dna_A + tseq("ghi")           // other possible variant
  });

  run_fwd_search_entry(tseq("ghi").size(), tseq("abcde").size(), tseq("bcde").size(),
                       f_pos + tseq("f").size(), tseq("ghij").size(),
                       tseq("bcde") + dna_A + tseq("ghij"));

  EXPECT_THAT(m_push_entries,
              ElementsAre(
                  // other possible variant
                  FwdPushSearchEntry(tseq("ghi").size(), f_pos + tseq("f").size(),
                                     dna_G + dna_A + tseq("ghij"), dna_G + dna_A + tseq("ghi"))));

  EXPECT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(tseq("a").size(), tseq("bcde") + dna_A + tseq("ghij"),
                                     f_pos + tseq("fghij").size())));
}

TEST_P(rejoin_test, rejoin_and_search_more2) {
  m_options.min_overlap = tseq("ghi").size();
  use_ref_parts({{0, tseq("abcde") + dna_T + tseq("fghijklmnopqrstuvwxyz")}});
  aoffset_t f_pos = (tseq("abcde") + dna_T).size();
  use_reads({
      tseq("bcde") + dna_A + tseq("ghij"),       // rejoining variant
      dna_T + tseq("cde") + dna_A + tseq("ghi")  // other possible variant
  });

  run_fwd_search_entry(tseq("ghi").size(), tseq("abcde").size(), tseq("bcde").size(),
                       f_pos + tseq("f").size(), tseq("ghij").size(),
                       tseq("bcde") + dna_A + tseq("ghij"));

  EXPECT_THAT(m_push_entries,
              ElementsAre(FwdPushSearchEntry(tseq("ghi").size(), f_pos + tseq("f").size(),
                                             dna_T + tseq("cde") + dna_A + tseq("ghij"),
                                             dna_T + tseq("cde") + dna_A + tseq("ghi"))));

  EXPECT_THAT(m_assemblies,
              ElementsAre(AssemblyIs(tseq("a").size(), tseq("bcde") + dna_A + tseq("ghij"),
                                     f_pos + tseq("fghij").size())));
}

INSTANTIATE_TEST_CASE_P(fwd_rejoin_test, rejoin_test, ::testing::Values(false /* not rev_comp */));
INSTANTIATE_TEST_CASE_P(rev_rejoin_test, rejoin_test, ::testing::Values(true /* rev_comp */));

}  // namespace discovery
}  // namespace variants

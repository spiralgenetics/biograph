#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/variants/discovery/discovery_testutil.h"
#include "modules/variants/discovery/state.h"

#include <fstream>

using namespace testing;
using namespace dna_testutil;

namespace variants {
namespace discovery {

class state_test : public discovery_test, public WithParamInterface<bool /* rev comp */> {
 public:
  state_test() { m_rev_comp = GetParam(); }

  void init_state() {
    init_discovery();
    add_ref_without_search();
    m_view = fwd_view();
  }

  void output_assembly(assembly_ptr a) {
    if (m_rev_comp) {
      reverse_assembly_in_place(a.get(), m_options.readmap, m_options.scaffold->end_pos());
    }
    m_st->output_assembly(std::move(a), true /* walk for additional variants */);
  }

  using ploids_remaining_t = std::vector<
      std::pair<std::pair<aoffset_t /* start */, aoffset_t /* end */>, int /* ploids left */>>;
  ploids_remaining_t ploids_remaining() {
    ploids_remaining_t result;
    for (const auto& elem : m_st->m_ploids_remaining) {
      if (m_rev_comp) {
        result.emplace_back(std::make_pair(m_view->reverse_offset(elem.first.upper()),
                                           m_view->reverse_offset(elem.first.lower())),
                            elem.second);
      } else {
        result.emplace_back(std::make_pair(elem.first.lower(), elem.first.upper()), elem.second);
      }
    }
    if (m_rev_comp) {
      std::reverse(result.begin(), result.end());
    }
    return result;
  }

  void output_assembly(aoffset_t left_offset, aoffset_t left_anchor_len, aoffset_t right_offset,
                       aoffset_t right_anchor_len, dna_slice seq) {
    assembly_ptr a = make_unique<assembly>();

    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    a->left_anchor_len = left_anchor_len;
    a->right_anchor_len = right_anchor_len;

    output_assembly(std::move(a));
  }

  view_t* m_view = nullptr;
};

MATCHER_P2(OffsetInfoIs, ploids_remaining, ref_remaining_limit, "") {
  if (arg.ploids_remaining != ploids_remaining) {
    return false;
  }
  if (aoffset_t(arg.ref_remaining_limit) != aoffset_t(ref_remaining_limit)) {
    return false;
  }
  return true;
}

Matcher<const offset_info&> OffsetInfoHasNoPloids() {
  return Field(&offset_info::ploids_remaining, Le(0));
}

TEST_P(state_test, ploids_remaining_init) {
  m_options.bidir_max_ploids = 123;
  use_ref_parts({{0, tseq("abcdefg")}});
  use_reads({tseq("abcdefg")});
  init_state();
  EXPECT_THAT(ploids_remaining(), ElementsAre(Pair(Pair(0, tseq("abcdefg").size()), 123)));

  EXPECT_THAT(m_view->get_offset_info(0, true /* fwd */),
              OffsetInfoIs(123, tseq("abcdefg").size()));

  EXPECT_THAT(m_view->get_offset_info(tseq("abcdefg").size(), false /* rev */),
              OffsetInfoIs(123, 0));
}

TEST_P(state_test, ploids_remaining_decrease) {
  m_options.bidir_max_ploids = 123;
  use_ref_parts({{0, tseq("abcdefg")}});
  use_reads({tseq("abcdefg")});
  init_state();
  output_assembly(tseq("a").size(), tseq("bc").size(), tseq("abcdef").size(), tseq("ef").size(),
                  tseq("bc") + dna_T + tseq("ef"));
  EXPECT_THAT(ploids_remaining(),
              ElementsAre(  //
                  Pair(Pair(0, tseq("ab").size() - 1), 123),
                  Pair(Pair(tseq("ab").size(), tseq("abcde").size()), 122),
                  Pair(Pair(tseq("abcde").size() + 1, tseq("abcdefg").size()), 123)));

  EXPECT_THAT(m_view->get_offset_info(0, true /* fwd */),
              OffsetInfoIs(123, tseq("abcdefg").size()));

  EXPECT_THAT(m_view->get_offset_info(tseq("abcdefg").size(), false /* rev */),
              OffsetInfoIs(123, 0));

  EXPECT_THAT(m_view->get_offset_info(tseq("ab").size(), true /* fwd */),
              OffsetInfoIs(122, tseq("abcdefg").size()));
  EXPECT_THAT(m_view->get_offset_info(tseq("abcde").size(), false /* rev */), OffsetInfoIs(122, 0));
}

TEST_P(state_test, use_all_ploids) {
  m_options.bidir_max_ploids = 1;
  use_ref_parts({{0, tseq("abcdefg")}});
  use_reads({tseq("abcdefg")});
  init_state();
  output_assembly(tseq("a").size(), tseq("bc").size(), tseq("abcdef").size(), tseq("ef").size(),
                  tseq("bc") + dna_T + tseq("ef"));
  EXPECT_THAT(ploids_remaining(),
              ElementsAre(  //
                  Pair(Pair(0, tseq("ab").size() - 1), 1),
                  Pair(Pair(tseq("abcde").size() + 1, tseq("abcdefg").size()), 1)));

  EXPECT_THAT(m_view->get_offset_info(0, true /* fwd */), OffsetInfoIs(1, tseq("ab").size() - 1));
  EXPECT_THAT(m_view->get_offset_info(tseq("ab").size() - 1, true /* fwd */),
              OffsetInfoIs(1, tseq("ab").size() - 1));
  EXPECT_THAT(m_view->get_offset_info(tseq("ab").size(), true /* fwd */), OffsetInfoHasNoPloids());

  EXPECT_THAT(m_view->get_offset_info(tseq("abcde").size(), false /* rev */),
              OffsetInfoHasNoPloids());
  EXPECT_THAT(m_view->get_offset_info(tseq("abcde").size() + 1, false /* rev */),
              OffsetInfoIs(1, tseq("abcde").size() + 1));
  EXPECT_THAT(m_view->get_offset_info(tseq("abcdefg").size(), false /* rev */),
              OffsetInfoIs(1, tseq("abcde").size() + 1));
}

TEST_P(state_test, use_some_ploids) {
  m_options.bidir_max_ploids = 2;
  use_ref_parts({{0, tseq("abcdefghi")}});
  use_reads({tseq("abcdefghi")});
  init_state();
  output_assembly(tseq("a").size(), tseq("bc").size(), tseq("abcdefgh").size(), tseq("gh").size(),
                  tseq("bc") + dna_T + tseq("gh"));
  output_assembly(tseq("ab").size(), tseq("cd").size(), tseq("abcdefg").size(), tseq("fg").size(),
                  tseq("cd") + dna_T + tseq("fg"));
  EXPECT_THAT(ploids_remaining(),
              ElementsAre(  //
                  Pair(Pair(0, tseq("ab").size() - 1), 2),
                  Pair(Pair(tseq("ab").size(), tseq("abc").size() - 1), 1),
                  Pair(Pair(tseq("abcdef").size() + 1, tseq("abcdefg").size()), 1),
                  Pair(Pair(tseq("abcdefg").size() + 1, tseq("abcdefghi").size()), 2)));

  EXPECT_THAT(m_view->get_offset_info(0, true /* fwd */), OffsetInfoIs(2, tseq("abc").size() - 1));
  EXPECT_THAT(m_view->get_offset_info(tseq("ab").size() - 1, true /* fwd */),
              OffsetInfoIs(2, tseq("abc").size() - 1));
  EXPECT_THAT(m_view->get_offset_info(tseq("ab").size(), true /* fwd */),
              OffsetInfoIs(1, tseq("abc").size() - 1));
  EXPECT_THAT(m_view->get_offset_info(tseq("abc").size() - 1, true /* fwd */),
              OffsetInfoIs(1, tseq("abc").size() - 1));
  EXPECT_THAT(m_view->get_offset_info(tseq("abc").size(), true /* fwd */), OffsetInfoHasNoPloids());

  EXPECT_THAT(m_view->get_offset_info(tseq("abcdefghi").size(), false /* rev */),
              OffsetInfoIs(2, tseq("abcdef").size() + 1));

  EXPECT_THAT(m_view->get_offset_info(tseq("abcdefg").size() + 1, false /* rev */),
              OffsetInfoIs(2, tseq("abcdef").size() + 1));

  EXPECT_THAT(m_view->get_offset_info(tseq("abcdefg").size(), false /* rev */),
              OffsetInfoIs(1, tseq("abcdef").size() + 1));

  EXPECT_THAT(m_view->get_offset_info(tseq("abcdef").size() + 1, false /* rev */),
              OffsetInfoIs(1, tseq("abcdef").size() + 1));

  EXPECT_THAT(m_view->get_offset_info(tseq("abcdef").size(), false /* rev */),
              OffsetInfoHasNoPloids());
}

INSTANTIATE_TEST_CASE_P(fwd_state_test, state_test, ::testing::Values(false /* not rev_comp */));
INSTANTIATE_TEST_CASE_P(rev_state_test, state_test, ::testing::Values(true /* rev_comp */));

}  // namespace discovery
}  // namespace variants

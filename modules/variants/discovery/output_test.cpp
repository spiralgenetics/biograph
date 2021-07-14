#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/variants/discovery/discovery.h"
#include "modules/variants/discovery/discovery_testutil.h"
#include "modules/variants/discovery/state.h"

#include <fstream>

using namespace testing;
using namespace dna_testutil;

namespace variants {
namespace discovery {

class output_test : public discovery_test {};

TEST_F(output_test, output_fwd) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd"), tseq("fghij")});
  init_discovery();

  std::shared_ptr<search_entry_t> e =
      std::make_shared<search_entry_t>(fwd_view(), 0, get_seqset_range(tseq("fghij")));
  e->left_anchor_len = tseq("abcd").size();
  e->seq = tseq("abcd") + dna_A;
  m_st->output_join_ref(e, tseq("abcde").size());

  ASSERT_THAT(m_assemblies, SizeIs(1)) << PrintToString(m_assemblies);
  const auto& a = *m_assemblies.begin();

  EXPECT_EQ(a.left_offset, 0);
  EXPECT_EQ(a.left_anchor_len, tseq("abcd").size());
  EXPECT_EQ(a.right_offset, tseq("abcdefghij").size());
  EXPECT_EQ(a.right_anchor_len, tseq("fghij").size());
  EXPECT_EQ(a.seq, tseq("abcd") + dna_A + tseq("fghij"));
}

TEST_F(output_test, output_rev) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd"), tseq("fghij")});
  init_discovery();

  std::shared_ptr<search_entry_t> e =
      std::make_shared<search_entry_t>(rev_view(), m_st->reverse_offset(tseq("abcdefghij").size()),
                                       get_seqset_range(tseq_rc("abcd")));
  e->left_anchor_len = tseq_rc("fghij").size();
  e->seq = tseq_rc("fghij") + dna_A.rev_comp();
  m_st->output_join_ref(e, m_st->reverse_offset(tseq("abcd").size()));

  ASSERT_THAT(m_assemblies, SizeIs(1)) << PrintToString(m_assemblies);
  const auto& a = *m_assemblies.begin();

  EXPECT_EQ(a.left_offset, 0);
  EXPECT_EQ(a.left_anchor_len, tseq("abcd").size());
  EXPECT_EQ(a.right_offset, tseq("abcdefghij").size());
  EXPECT_EQ(a.right_anchor_len, tseq("fghij").size());
  EXPECT_EQ(a.seq, tseq("abcd") + dna_A + tseq("fghij"));
}

}  // namespace discovery
}  // namespace variants

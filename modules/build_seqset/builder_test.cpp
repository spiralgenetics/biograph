#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/seqset_testutil.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace testing;
using namespace dna_testutil;

}  // namespace

class builder_test : public TestWithParam<unsigned /* partition depth */> {
 public:
  void SetUp() override { g_seqset_build_partition_depth = GetParam(); }

  void verify_seqset(const seqset& ss) {
    dna_sequence last;
    for (size_t i = 0; i < ss.size(); ++i) {
      auto e = ss.ctx_entry(i);
      auto seq = e.sequence();
      EXPECT_GT(seq.size(), 0);
      EXPECT_EQ(int(last.shared_prefix_length(seq)), int(ss.entry_shared(i)))
          << "Last:\n"
          << last << "\nSeq:\n"
          << seq;

      if (last.size() != 0) {
        EXPECT_EQ(last.compare_to(seq), dna_compare_result::FIRST_IS_LESS)
            << "Last:\n"
            << last << "\nSeq:\n"
            << seq;
      }

      auto popped = ss.entry_pop_front(i);
      EXPECT_GE(ss.entry_size(popped), ss.entry_size(i) - 1) << seq;
      last = seq;
    }
  }

  std::set<dna_sequence> seqset_entries(const seqset& ss) {
    std::set<dna_sequence> result;

    for (size_t i = 0; i < ss.size(); ++i) {
      EXPECT_TRUE(result.insert(ss.ctx_entry(i).sequence()).second);
    }
    return result;
  }
};

TEST_P(builder_test, seq1) {
  auto ss_f = seqset_for_reads({tseq("a")});
  const auto& ss = ss_f->get_seqset();

  verify_seqset(ss);
  EXPECT_THAT(
      seqset_entries(ss),
      ElementsAreArray({dna_sequence("AAAATTAC"), dna_sequence("AAATTAC"),
                        dna_sequence("AATTAC"), dna_sequence("AATTTTAG"),
                        dna_sequence("AC"), dna_sequence("AG"),
                        dna_sequence("ATTAC"), dna_sequence("ATTTTAG"),
                        dna_sequence("CTAAAATTAC"), dna_sequence("GTAATTTTAG"),
                        dna_sequence("TAAAATTAC"), dna_sequence("TAATTTTAG"),
                        dna_sequence("TAC"), dna_sequence("TAG"),
                        dna_sequence("TTAC"), dna_sequence("TTAG"),
                        dna_sequence("TTTAG"), dna_sequence("TTTTAG")}));
}

TEST_P(builder_test, seq2) {
  auto ss_f = seqset_for_reads({tseq("abcdefg")});
  const auto& ss = ss_f->get_seqset();

  verify_seqset(ss);
  EXPECT_EQ(129, ss.size());
}

TEST_P(builder_test, seq3) {
  auto ss_f = seqset_for_reads({tseq("abcd"), tseq("cdef"), tseq_rc("efgh")});
  const auto& ss = ss_f->get_seqset();

  verify_seqset(ss);
  EXPECT_EQ(152, ss.size());
}

TEST_P(builder_test, seq4) {
  auto ss_f =
      seqset_for_reads({tseq("ab"), tseq("bc"), tseq("cd"), tseq("be")});
  const auto& ss = ss_f->get_seqset();

  verify_seqset(ss);
  EXPECT_EQ(91, ss.size());
}
TEST_P(builder_test, seq5) {
  auto ss_f =
      seqset_for_reads({tseq("AB"), tseq("BC"), tseq("CD"), tseq("BE")});
  const auto& ss = ss_f->get_seqset();

  verify_seqset(ss);
  EXPECT_EQ(91, ss.size());
}

TEST_P(builder_test, seq6) {
  auto ss_f = seqset_for_reads({tseq("abc"), tseq("cde")});
  const auto& ss = ss_f->get_seqset();
  verify_seqset(ss);
  EXPECT_TRUE(ss.find(tseq("abc")).valid());
  EXPECT_TRUE(ss.find(tseq_rc("abc")).valid());
  EXPECT_TRUE(ss.find(tseq_rc("cde")).valid());
  EXPECT_TRUE(ss.find(tseq("cde")).valid());
  EXPECT_EQ(89, ss.size());
}

TEST_P(builder_test, seq7) {
  auto ss_f = seqset_for_reads({tseq("abc"), tseq("efg")});
  const auto& ss = ss_f->get_seqset();
  verify_seqset(ss);
  EXPECT_TRUE(ss.find(tseq("abc")).valid());
  EXPECT_TRUE(ss.find(tseq_rc("abc")).valid());
  EXPECT_TRUE(ss.find(tseq_rc("efg")).valid());
  EXPECT_TRUE(ss.find(tseq("efg")).valid());
  EXPECT_EQ(99, ss.size());
}

INSTANTIATE_TEST_CASE_P(  //
    builder_test, builder_test, ::testing::Values(1, 2, 3, 4));

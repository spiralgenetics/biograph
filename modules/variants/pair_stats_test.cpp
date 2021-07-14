#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/variants/pair_stats.h"
#include "modules/variants/ref_map.h"

using namespace testing;
using namespace dna_testutil;
using namespace variants;

class pair_stats_test : public Test {
 protected:
  pair_stats_test() {
    if (!m_ref) {
      m_ref = create_reference({tseq("abcdefghijklmnopqrstuvwxyz0123456789"),
                                tseq("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")});
    }
  }

  void use_reads(
      const std::vector<std::pair<dna_sequence, dna_sequence>>& pairs) {
    std::vector<dna_sequence> all_reads;
    for (const auto& pair : pairs) {
      all_reads.push_back(pair.first);
      all_reads.push_back(pair.second);
    }
    m_seqset = seqset_for_reads(all_reads);
    m_readmap = readmap_for_reads(m_seqset, pairs, {});

    m_rmap.emplace(&*m_seqset, m_ref.get());
    m_rmap->build();

    m_pair_stats.emplace(m_seqset.get(), m_readmap.get(), m_ref.get(), &m_rmap.get());
    m_pair_stats->calc_stats();
  }

  boost::optional<ref_map> m_rmap;
  boost::optional<pair_stats> m_pair_stats;
  static std::unique_ptr<reference> m_ref;
  std::shared_ptr<seqset> m_seqset;
  std::unique_ptr<readmap> m_readmap;
};

std::unique_ptr<reference> pair_stats_test::m_ref;

TEST_F(pair_stats_test, empty) {
  use_reads({{tseq("abcde"), tseq("ABCDE")}});
  EXPECT_FALSE(m_pair_stats->found_pairs());
}

TEST_F(pair_stats_test, simple) {
  use_reads({{tseq("abcde"), tseq_rc("hijkl")}});
  ASSERT_TRUE(m_pair_stats->found_pairs());
  EXPECT_EQ(tseq("abcdefghijkl").size(), m_pair_stats->median_pair_offset());
}

TEST_F(pair_stats_test, simple_rc) {
  use_reads({{tseq_rc("hijkl"), tseq("abcde")}});
  ASSERT_TRUE(m_pair_stats->found_pairs());
  EXPECT_EQ(tseq("abcdefghijkl").size(), m_pair_stats->median_pair_offset());
}

TEST_F(pair_stats_test, multi) {
  use_reads({
      {tseq("abcde"), tseq_rc("hijkl")},  // distance 12
      {tseq("bcdef"), tseq_rc("jklmn")},  // distance 13
      {tseq("cdefg"), tseq_rc("lmnop")}  // distance 14
  });
  ASSERT_TRUE(m_pair_stats->found_pairs());
  EXPECT_GT(m_pair_stats->median_pair_offset(), tseq("abcdefghijkl").size());
  EXPECT_LT(m_pair_stats->median_pair_offset(),
            tseq("cdefghijklmnop").size());
}

TEST_F(pair_stats_test, reverse) {
  use_reads({{tseq_rc("abcde"), tseq("hijkl")}});
  EXPECT_EQ(-int64_t(tseq("fg").size()), m_pair_stats->median_pair_offset());
}

TEST_F(pair_stats_test, reverse_rc) {
  use_reads({{tseq("hijkl"), tseq_rc("abcde")}});
  ASSERT_TRUE(m_pair_stats->found_pairs());
  EXPECT_EQ(-int64_t(tseq("fg").size()), m_pair_stats->median_pair_offset());
}

TEST_F(pair_stats_test, disregard_bad_direction) {
  use_reads({{tseq("abcde"), tseq("hijkl")}});
  EXPECT_FALSE(m_pair_stats->found_pairs());
}

TEST_F(pair_stats_test, disregard_bad_direction_rc) {
  use_reads({{tseq_rc("abcde"), tseq_rc("hijkl")}});
  EXPECT_FALSE(m_pair_stats->found_pairs());
}

TEST_F(pair_stats_test, disregard_ambiguous) {
  use_reads({{tseq_rc("abcde"), tseq("012345")}});
  EXPECT_FALSE(m_pair_stats->found_pairs());
}

TEST_F(pair_stats_test, disregard_ambiguous_rc) {
  use_reads({{tseq("abcde"), tseq_rc("012345")}});
  EXPECT_FALSE(m_pair_stats->found_pairs());
}

TEST_F(pair_stats_test, disregard_cross_extent) {
  use_reads({{tseq("abcde"), tseq_rc("HIJKL")}});
  EXPECT_FALSE(m_pair_stats->found_pairs());
}

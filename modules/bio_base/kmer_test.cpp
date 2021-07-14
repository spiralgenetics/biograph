#include "modules/bio_base/kmer.h"

#include <gtest/gtest.h>

using namespace testing;

class kmer_test
    : public TestWithParam<
          bool /* true = use kmer_str_view, false = use kmer_view */> {
 public:
  kmer_test() { m_use_kmer_str = GetParam(); }

  bool m_use_kmer_str;
};

TEST_P(kmer_test, iterator) {
  const size_t kmer_size = 9;
  const std::string seqstr = "ACGTACGTACGTACGTACGT";
  SPLOG("%s", seqstr.c_str());
  dna_sequence seq(seqstr);

  std::vector<std::string> kmer_strs;
  std::vector<kmer_t> kmers;

  if (m_use_kmer_str) {
    for (const auto& kmer : kmer_str_view(seqstr, kmer_size)) {
      dna_sequence kmer_seq(kmer, kmer_size);
      auto kmer_str = kmer_seq.as_string();
      SPLOG("%s: 0x%016zx", kmer_str.c_str(), kmer);
      kmer_strs.push_back(kmer_str);
      kmers.push_back(kmer);
    }
  } else {
    for (const auto& kmer : kmer_view(seq, kmer_size)) {
      dna_sequence kmer_seq(kmer, kmer_size);
      auto kmer_str = kmer_seq.as_string();
      SPLOG("%s: 0x%016zx", kmer_str.c_str(), kmer);
      kmer_strs.push_back(kmer_str);
      kmers.push_back(kmer);
    }
  }

  EXPECT_EQ(seq.size() - kmer_size + 1, kmer_strs.size());
  EXPECT_EQ("ACGTACGTA", kmer_strs[0]);
  EXPECT_EQ("CGTACGTAC", kmer_strs[1]);
  EXPECT_EQ("GTACGTACG", kmer_strs[2]);
  EXPECT_EQ("TACGTACGT", kmer_strs[3]);
  EXPECT_EQ("ACGTACGTA", kmer_strs[4]);
  EXPECT_EQ("CGTACGTAC", kmer_strs[5]);
  EXPECT_EQ("GTACGTACG", kmer_strs[6]);
  EXPECT_EQ("TACGTACGT", kmer_strs[7]);
  EXPECT_EQ("ACGTACGTA", kmer_strs[8]);
  EXPECT_EQ("CGTACGTAC", kmer_strs[9]);
  EXPECT_EQ("GTACGTACG", kmer_strs[10]);
  EXPECT_EQ("TACGTACGT", kmer_strs[11]);

  EXPECT_EQ(0x0000000000006c6c, kmers[0]);
  EXPECT_EQ(0x000000000001b1b1, kmers[1]);
  EXPECT_EQ(0x000000000002c6c6, kmers[2]);
  EXPECT_EQ(0x0000000000031b1b, kmers[3]);
  EXPECT_EQ(0x0000000000006c6c, kmers[4]);
  EXPECT_EQ(0x000000000001b1b1, kmers[5]);
  EXPECT_EQ(0x000000000002c6c6, kmers[6]);
  EXPECT_EQ(0x0000000000031b1b, kmers[7]);
  EXPECT_EQ(0x0000000000006c6c, kmers[8]);
  EXPECT_EQ(0x000000000001b1b1, kmers[9]);
  EXPECT_EQ(0x000000000002c6c6, kmers[10]);
  EXPECT_EQ(0x0000000000031b1b, kmers[11]);
}

TEST_P(kmer_test, iterator_long) {
  const size_t kmer_size = 17;
  const std::string seqstr = "ACGTACGTACGTACGTACGT";
  SPLOG("%s", seqstr.c_str());
  dna_sequence seq(seqstr);

  std::vector<std::string> kmer_strs;
  std::vector<kmer_t> kmers;

  if (m_use_kmer_str) {
    for (const auto& kmer : kmer_str_view(seqstr, kmer_size)) {
      dna_sequence kmer_seq(kmer, kmer_size);
      auto kmer_str = kmer_seq.as_string();
      SPLOG("%s: 0x%016zx", kmer_str.c_str(), kmer);
      kmer_strs.push_back(kmer_str);
      kmers.push_back(kmer);
    }
  } else {
    for (const auto& kmer : kmer_view(seq, kmer_size)) {
      dna_sequence kmer_seq(kmer, kmer_size);
      auto kmer_str = kmer_seq.as_string();
      SPLOG("%s: 0x%016zx", kmer_str.c_str(), kmer);
      kmer_strs.push_back(kmer_str);
      kmers.push_back(kmer);
    }
  }

  EXPECT_EQ(seq.size() - kmer_size + 1, kmer_strs.size());
  EXPECT_EQ("ACGTACGTACGTACGTA", kmer_strs[0]);
  EXPECT_EQ("CGTACGTACGTACGTAC", kmer_strs[1]);
  EXPECT_EQ("GTACGTACGTACGTACG", kmer_strs[2]);
  EXPECT_EQ("TACGTACGTACGTACGT", kmer_strs[3]);

  EXPECT_EQ(0x000000006c6c6c6c, kmers[0]);
  EXPECT_EQ(0x00000001b1b1b1b1, kmers[1]);
  EXPECT_EQ(0x00000002c6c6c6c6, kmers[2]);
  EXPECT_EQ(0x000000031b1b1b1b, kmers[3]);
}

INSTANTIATE_TEST_CASE_P(kmer_view_test, kmer_test,
                        ::testing::Values(false /* kmer_view */));
INSTANTIATE_TEST_CASE_P(kmer_str_view_test, kmer_test,
                        ::testing::Values(true /* kmer_str_view */));

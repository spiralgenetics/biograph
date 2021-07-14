#include "modules/bio_base/dna_testutil.h"
#include "modules/build_seqset/kmer_counter.h"
#include "modules/io/parallel.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <random>

namespace build_seqset {

using namespace testing;
using namespace dna_testutil;

// Make sure k_unused_entry can never be encountered in practice.
TEST(kmer_count_table_test, unused_entry) {
  kmer_t k = kmer_count_table<uint8_t>::k_unused_entry & kmer_count_table<uint8_t>::k_kmer_mask;
  bool flipped;

  kmer_t canon = canonicalize(k, 30, flipped);
  EXPECT_NE(canon, k);
  EXPECT_TRUE(flipped);
}

struct expected_data {
  uint32_t fwd_count = 0;
  uint32_t rev_count = 0;
  bool fwd_flag = false;
  bool rev_flag = false;
};

TEST(kmer_count_table_test, flags) {
  using table_type = kmer_count_table<uint32_t>;
  constexpr unsigned k_kmer_size = 30;

  std::random_device rand_dev;
  std::mt19937_64 rand_source(rand_dev());
  tracked_vector<kmer_t> input_kmers(track_alloc("kmer_count_table_test:input_kmers"));
  std::map<kmer_t, expected_data> expected;

  table_type table(400, "kmer_count_table_test");

  for (unsigned i = 0; i < 200; ++i) {
    kmer_t rand_kmer = rand_source() & ((1ULL << (k_kmer_size * 2)) - 1);
    kmer_t canon = canonicalize(rand_kmer, k_kmer_size);
    expected_data& e = expected[canon];

    unsigned nentries = std::uniform_int_distribution<unsigned>(0, 100000)(rand_source);
    unsigned nflipped = std::uniform_int_distribution<unsigned>(0, nentries)(rand_source);

    for (unsigned j = 0; j < nentries; ++j) {
      bool flipped = j < nflipped;
      kmer_t k = flipped ? rev_comp(canon, k_kmer_size) : canon;
      uint32_t& count = flipped ? e.rev_count : e.fwd_count;
      bool& fwd_flag = flipped ? e.rev_flag : e.fwd_flag;
      bool& rev_flag = flipped ? e.fwd_flag : e.rev_flag;

      count++;

      if (std::uniform_int_distribution<>(0, 1000)(rand_source) == 0) {
        if (std::uniform_int_distribution<>(0, 10)(rand_source) == 0) {
          k |= table_type::k_fwd_flag;
          fwd_flag = true;
        }
        if (std::uniform_int_distribution<>(0, 10)(rand_source) == 0) {
          k |= table_type::k_rev_flag;
          rev_flag = true;
        }
      }
      input_kmers.push_back(k);
    }
  }

  std::shuffle(input_kmers.begin(), input_kmers.end(), rand_source);

  parallel_for(0, input_kmers.size(), [&input_kmers, &table](size_t start, size_t limit) {
    auto end = input_kmers.begin() + limit;
    for (auto it = input_kmers.begin() + start; it != end; ++it) {
      kmer_t k_and_flags = *it;
      kmer_t k = k_and_flags & table_type::k_kmer_mask;
      bool flipped;
      kmer_t canon = canonicalize(k, k_kmer_size, flipped);

      table.increment(canon, flipped, k_and_flags & table_type::k_fwd_flag,
                      k_and_flags & table_type::k_rev_flag);
    }
  });

  for (const auto& exp : expected) {
    kmer_t kmer = exp.first;
    const auto& elem = table.get(kmer);
    EXPECT_TRUE(elem) << "Couldn't find kmer " << printstring("%lX", kmer);
    if (!elem) {
      continue;
    }
    EXPECT_EQ(kmer, elem.kmer());
    const expected_data& e = exp.second;
    EXPECT_EQ(e.fwd_count, elem.fwd_count);
    EXPECT_EQ(e.rev_count, elem.rev_count);
    EXPECT_EQ(e.fwd_flag, elem.fwd_flag());
    EXPECT_EQ(e.rev_flag, elem.rev_flag());
  }

  table.compact();
  for (auto it = table.begin(); it != table.end(); ++it) {
    ASSERT_TRUE(*it);
    auto e_it = expected.find(it->kmer());
    EXPECT_TRUE(e_it != expected.end()) << "Couldn't find kmer " << printstring("%lX", it->kmer());
    if (e_it == expected.end()) {
      continue;
    }

    const auto& elem = *it;
    const expected_data& expect = e_it->second;

    EXPECT_EQ(expect.fwd_count, elem.fwd_count);
    EXPECT_EQ(expect.rev_count, elem.rev_count);
    EXPECT_EQ(expect.fwd_flag, elem.fwd_flag());
    EXPECT_EQ(expect.rev_flag, elem.rev_flag());

    expected.erase(e_it);
  }

  EXPECT_THAT(expected, IsEmpty()) << "Some results missing from count";
}

}  // namespace build_seqset

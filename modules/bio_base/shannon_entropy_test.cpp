#include "modules/bio_base/shannon_entropy.h"
#include "modules/bio_base/dna_testutil.h"

#include <gtest/gtest.h>

TEST(shannon_entropy_test, random_gen) {
  std::random_device true_rand;
  std::mt19937 rand_source(true_rand());

  constexpr unsigned k_seq_len = 200;
  constexpr unsigned k_num_seqs = 100;
  constexpr unsigned k_target_len = 70;

  uint64_t entropy_sum = 0;
  uint64_t length_sum = 0;
  uint64_t entropy_sample_count = 0;

  for (unsigned i = 0; i < k_num_seqs; ++i) {
    dna_sequence seq = rand_dna_sequence(rand_source, k_seq_len);

    shannon_entropy e(k_target_len);

    unsigned idx = 0;
    unsigned last_needed = 0;
    for (dna_base b : seq) {
      e.push_front(b);
      auto needed = e.length_needed();

      auto entropy = e.calc_entropy();

      if (needed) {
        EXPECT_GT(idx, k_target_len);
      } else {
        // True random sources will occasionally supply low-entropy
        // sequences, so be pretty tolerant.
        EXPECT_LT(idx, k_target_len * 2);
      }

      idx++;

      if (!needed) {
        continue;
      }

      if (last_needed) {
        EXPECT_LE(*needed, (last_needed + 1));
      }

      last_needed = *needed;
      entropy_sum += entropy;
      length_sum += *needed;
      entropy_sample_count++;
    }
  }

  double entropy_avg = entropy_sum * 1. / entropy_sample_count;
  EXPECT_LT(entropy_avg, k_target_len);
  EXPECT_GT(entropy_avg, k_target_len * 0.95);

  double length_avg = length_sum * 1. / entropy_sample_count;
  EXPECT_GT(length_avg, k_target_len);
  EXPECT_LT(length_avg, k_target_len * 1.2);
}

TEST(shannon_entropy_test, simple) {
  dna_sequence seq("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

  shannon_entropy e(5);
  e.push_front(seq);

  auto needed = e.length_needed();
  EXPECT_FALSE(needed);

  EXPECT_EQ(0, e.calc_entropy());
}

TEST(shannon_entropy_test, repetitive) {
  dna_sequence seq(
      "ACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACT"
      "ACTACTACTACTACTACTACTACTACTA");

  shannon_entropy e(100);
  e.push_front(seq);

  auto needed = e.length_needed();
  EXPECT_FALSE(needed);

  EXPECT_EQ(25, e.calc_entropy());
}

TEST(shannon_entropy_test, rand_then_simple) {
  // Initially, fill it with random with high entopy.
  std::random_device true_rand;
  std::mt19937 rand_source(true_rand());
  dna_sequence seq(rand_dna_sequence(rand_source, 300));

  shannon_entropy e(150);
  e.push_front(seq);

  auto needed = e.length_needed();
  ASSERT_TRUE(needed);
  EXPECT_GT(*needed, 150);
  EXPECT_LT(*needed, 150 * 1.2);
  unsigned max_pushed_needed = *needed;

  auto entropy = e.calc_entropy();
  EXPECT_GT(entropy, 150 * 0.9);
  EXPECT_LT(entropy, 150);

  // Allow entropy to potentially get slightly bigger when start
  // pushing A's, but afterwards should always decrease.
  auto last_entropy = entropy;
  // Next, fill it with repetetive As. This should make the entropy go
  // down.
  for (unsigned i = 0; i < 300; ++i) {
    e.push_front(dna_base('A'));
    auto new_entropy = e.calc_entropy();
    needed = e.length_needed();
    if (new_entropy < (entropy - 3)) {
      EXPECT_LE(new_entropy, last_entropy);
      last_entropy = new_entropy;

      if (needed) {
        EXPECT_GE(*needed, max_pushed_needed);
      }
    }
    if (needed) {
      max_pushed_needed = *needed;
    }
  }

  EXPECT_EQ(255, max_pushed_needed);
  EXPECT_FALSE(needed);

  entropy = e.calc_entropy();
  EXPECT_EQ(0, entropy);
}

#include "modules/variants/align.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

struct match_test_param {
  int left_seq_size;
  int left_scaffold_size;
  int match_size;
  int right_seq_size;
  int right_scaffold_size;

  friend std::ostream& operator<<(std::ostream& os,
                                  const match_test_param& param) {
    os << "seq: " << param.left_seq_size << " + " << param.match_size << " + "
       << param.right_seq_size << ", "
       << "scaffold: " << param.left_scaffold_size << " + " << param.match_size
       << " + " << param.right_scaffold_size;
    return os;
  }
};

dna_sequence repeat_dna(dna_sequence b, int len) {
  dna_sequence seq;
  for (int i = 0; i < len; ++i) {
    seq += b;
  }
  return seq;
}

class find_match_test : public TestWithParam<match_test_param> {
 public:
  find_match_test() {
    m_left_seq_size = GetParam().left_seq_size;
    m_left_scaffold_size = GetParam().left_scaffold_size;
    m_match_size = GetParam().match_size;
    m_right_seq_size = GetParam().right_seq_size;
    m_right_scaffold_size = GetParam().right_scaffold_size;
  }

 protected:
  int m_left_seq_size;
  int m_left_scaffold_size;
  int m_match_size;
  int m_right_seq_size;
  int m_right_scaffold_size;
};

TEST_P(find_match_test, matches_full) {
  dna_sequence seq(repeat_dna(dna_T, m_left_seq_size) +
                   repeat_dna(dna_C, m_match_size) +
                   repeat_dna(dna_T, m_right_seq_size));
  scaffold s(repeat_dna(dna_A, m_left_scaffold_size) +
             repeat_dna(dna_C, m_match_size) +
             repeat_dna(dna_A, m_right_scaffold_size));

  aoffset_t seq_match_start, scaffold_match_start;
  bool did_match = aligner::find_match(
      seq, s, m_match_size, &seq_match_start, &scaffold_match_start,
      aligner::anchor_type_t::ANCHORED_TO_BOTH);
  ASSERT_TRUE(did_match);

  EXPECT_EQ(seq_match_start, m_left_seq_size);
  EXPECT_EQ(scaffold_match_start, m_left_scaffold_size);
}

TEST_P(find_match_test, is_biggest_match) {
  dna_sequence seq(repeat_dna(dna_T, m_left_seq_size) +
                   repeat_dna(dna_C, m_match_size) +
                   repeat_dna(dna_T, m_right_seq_size));
  scaffold s(repeat_dna(dna_A, m_left_scaffold_size) +
             repeat_dna(dna_C, m_match_size) +
             repeat_dna(dna_A, m_right_scaffold_size));

  int match_len;
  aoffset_t seq_match_start, scaffold_match_start;
  assemble_options opts;
  opts.max_ref_align_bases = 1;
  int min_match_len;
  bool did_match = aligner::find_biggest_match(
      opts, seq, s, &match_len, &seq_match_start, &scaffold_match_start,
      &min_match_len, aligner::anchor_type_t::ANCHORED_TO_BOTH);
  ASSERT_TRUE(did_match);

  EXPECT_EQ(seq_match_start, m_left_seq_size);
  EXPECT_EQ(scaffold_match_start, m_left_scaffold_size);
  EXPECT_EQ(match_len, m_match_size);
}

TEST_P(find_match_test, no_biggest_match) {
  dna_sequence seq(repeat_dna(dna_T, m_left_seq_size) +
                   repeat_dna(dna_G, m_match_size) +
                   repeat_dna(dna_T, m_right_seq_size));
  scaffold s(repeat_dna(dna_A, m_left_scaffold_size) +
             repeat_dna(dna_C, m_match_size) +
             repeat_dna(dna_A, m_right_scaffold_size));

  int match_len;
  aoffset_t seq_match_start, scaffold_match_start;
  assemble_options opts;
  opts.max_ref_align_bases = 1;
  int min_match_len;
  bool did_match = aligner::find_biggest_match(
      opts, seq, s, &match_len, &seq_match_start, &scaffold_match_start,
      &min_match_len, aligner::anchor_type_t::ANCHORED_TO_BOTH);
  ASSERT_FALSE(did_match);
}

TEST_P(find_match_test, matches_sparse) {
  dna_sequence seq(repeat_dna(dna_T, m_left_seq_size) +
                   repeat_dna(dna_C, m_match_size) +
                   repeat_dna(dna_T, m_right_seq_size));
  scaffold s;
  s.add(m_left_scaffold_size, repeat_dna(dna_C, m_match_size));
  s.set_end_pos(m_left_scaffold_size + m_match_size + m_right_scaffold_size);

  aoffset_t seq_match_start, scaffold_match_start;
  bool did_match = aligner::find_match(
      seq, s, m_match_size, &seq_match_start, &scaffold_match_start,
      aligner::anchor_type_t::ANCHORED_TO_BOTH);
  ASSERT_TRUE(did_match);

  EXPECT_EQ(seq_match_start, m_left_seq_size);
  EXPECT_EQ(scaffold_match_start, m_left_scaffold_size);
}

TEST_P(find_match_test, doesnt_match_enough) {
  dna_sequence seq(repeat_dna(dna_T, m_left_seq_size) +
                   repeat_dna(dna_C, m_match_size) +
                   repeat_dna(dna_T, m_right_seq_size));
  scaffold s(repeat_dna(dna_A, m_left_scaffold_size) +
             repeat_dna(dna_C, m_match_size) +
             repeat_dna(dna_A, m_right_scaffold_size));

  aoffset_t seq_match_start, scaffold_match_start;
  bool did_match = aligner::find_match(
      seq, s, m_match_size + 1, &seq_match_start, &scaffold_match_start,
      aligner::anchor_type_t::ANCHORED_TO_BOTH);
  EXPECT_FALSE(did_match);
}

TEST_P(find_match_test, sparse_left_trunc) {
  dna_sequence seq(repeat_dna(dna_T, m_left_seq_size) +
                   repeat_dna(dna_C, m_match_size) +
                   repeat_dna(dna_T, m_right_seq_size));
  scaffold s;
  s.add(m_left_scaffold_size + 1, repeat_dna(dna_C, m_match_size - 1));
  s.set_end_pos(m_left_scaffold_size + m_match_size + m_right_scaffold_size);

  aoffset_t seq_match_start, scaffold_match_start;
  bool did_match = aligner::find_match(
      seq, s, m_match_size, &seq_match_start, &scaffold_match_start,
      aligner::anchor_type_t::ANCHORED_TO_BOTH);
  EXPECT_FALSE(did_match);
}

TEST_P(find_match_test, sparse_right_trunc) {
  dna_sequence seq(repeat_dna(dna_T, m_left_seq_size) +
                   repeat_dna(dna_C, m_match_size) +
                   repeat_dna(dna_T, m_right_seq_size));
  scaffold s;
  s.add(m_left_scaffold_size, repeat_dna(dna_C, m_match_size - 1));
  s.set_end_pos(m_left_scaffold_size + m_match_size + m_right_scaffold_size);

  aoffset_t seq_match_start, scaffold_match_start;
  bool did_match = aligner::find_match(
      seq, s, m_match_size, &seq_match_start, &scaffold_match_start,
      aligner::anchor_type_t::ANCHORED_TO_BOTH);
  EXPECT_FALSE(did_match);
}

std::vector<match_test_param> test_params() {
  std::vector<match_test_param> result;
  for (int left_seq_size : {0, 1, 100}) {
    for (int left_scaffold_size : {0, 1, 100}) {
      for (int right_seq_size : {0, 1, 100}) {
        for (int right_scaffold_size : {0, 1, 100}) {
          for (int match_size : {1, 29, 30, 31, 100}) {
            match_test_param p;
            p.left_seq_size = left_seq_size;
            p.left_scaffold_size = left_scaffold_size;
            p.right_seq_size = right_seq_size;
            p.right_scaffold_size = right_scaffold_size;
            p.match_size = match_size;
            result.push_back(p);
          }
        }
      }
    }
  }
  return result;
}

INSTANTIATE_TEST_CASE_P(find_match_tests, find_match_test,
                        ::testing::ValuesIn(test_params()));

}  // namespace variants

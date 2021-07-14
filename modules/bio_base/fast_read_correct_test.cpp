#include "modules/bio_base/fast_read_correct.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/io/parallel.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <set>
#include <unordered_set>

using namespace testing;
using namespace dna_testutil;

enum class error_behavior {
  // Replace bases with 'N' to indicate a read error.
  REPLACE_WITH_N,
  // Replace bases with different bases when injecting errors
  DIFFERENT_BASE
};

class fast_read_correct_test
    : public TestWithParam<std::tuple<int /* sequence size */, error_behavior>> {
 protected:
  void SetUp() override {
    std::tie(m_seq_size, m_read_error) = GetParam();
    dna_sequence long_test_seq = tseq("aBcDeFgHiJkLmNoPqRsTuVwXyZ");
    CHECK_GE(long_test_seq.size(), m_seq_size);
    dna_sequence seq = long_test_seq.subseq(0, m_seq_size);
    m_seq = seq.as_string();
    m_params.kmer_lookup_f = [this](kmer_t kmer, frc_kmer* kmer_info) {
      auto it = m_kmers.find(kmer);
      if (it != m_kmers.end()) {
        kmer_info->index = std::distance(m_kmers.begin(), it);
        kmer_info->flipped = false;
        return true;
      }
      it = m_kmers.find(rev_comp(kmer, m_params.kmer_size));
      if (it != m_kmers.end()) {
        kmer_info->index = std::distance(m_kmers.begin(), it);
        kmer_info->flipped = true;
        return true;
      }

      return false;
    };
  };

  void add_kmers() {
    CHECK_GE(m_seq.size(), m_params.kmer_size);
    for (size_t i = 0; i + m_params.kmer_size <= m_seq.size(); ++i) {
      std::string kmer = m_seq.substr(i, m_params.kmer_size);
      if (kmer.find('N') == std::string::npos) {
        dna_sequence kmer_seq(kmer);
        m_kmers.insert(make_kmer(kmer_seq.begin(), m_params.kmer_size));
      }
    }
  }

  std::string add_error(const std::string& orig, unsigned seq_offset) {
    CHECK_LT(seq_offset, m_seq.size());
    char error_base = m_seq[seq_offset];

    switch (m_read_error) {
      case error_behavior::REPLACE_WITH_N:
        error_base = 'N';
        break;
      case error_behavior::DIFFERENT_BASE:
        error_base = char(dna_base(int(dna_base(error_base)) ^ 1));
        break;
    }
    CHECK_GT(orig.size(), seq_offset);
    std::string err_seq = orig;
    err_seq[seq_offset] = error_base;
    return err_seq;
  }

  frc_output run_fast_read_correct(string_view input) {
    frc_output out = fast_read_correct(input, m_params);
    auto expected = expected_kmers(out.corrected);
    if (expected != out.kmers) {
      EXPECT_THAT(out.kmers, ContainerEq(expected_kmers(out.corrected)));
    }
    return out;
  }

  std::vector<frc_kmer> expected_kmers(dna_slice corrected) {
    std::vector<frc_kmer> result;
    for (kmer_t kmer : kmer_view(corrected, m_params.kmer_size)) {
      frc_kmer kmer_info;
      bool lookup_succeeded = m_params.kmer_lookup_f(kmer, &kmer_info);
      if (lookup_succeeded) {
        result.push_back(kmer_info);
      } else {
        EXPECT_TRUE(lookup_succeeded)
            << "Corrected: " << corrected << " kmer: " << kmer_str(kmer, m_params.kmer_size);
      }
    }

    return result;
  }

  error_behavior m_read_error;
  unsigned m_seq_size;
  std::string m_seq;
  std::unordered_set<kmer_t> m_kmers;
  frc_params m_params;
};

TEST_P(fast_read_correct_test, no_errors) {
  SCOPED_TRACE(m_seq.size());
  add_kmers();
  frc_output result = run_fast_read_correct(m_seq);
  EXPECT_EQ(result.corrected, m_seq);
  EXPECT_EQ(result.corrections, 0);
}

TEST_P(fast_read_correct_test, single) {
  SCOPED_TRACE(m_seq.size());
  add_kmers();
  for (unsigned i = 0; i < m_seq.size() - m_params.min_good_run; ++i) {
    std::string err_seq = add_error(m_seq, i);
    dna_sequence expected = m_seq;
    unsigned expected_corrections = 1;

    if (i < m_params.kmer_size && i >= m_seq.size() - m_params.kmer_size) {
      expected = dna_sequence();
      expected_corrections = 0;
    }

    SCOPED_TRACE(PrintToString(i) + PrintToString(err_seq));
    frc_output result = run_fast_read_correct(err_seq);
    EXPECT_EQ(result.corrected, expected) << " i: " << i;
    EXPECT_EQ(result.corrections, expected_corrections) << " i: " << i;
  }
}

TEST_P(fast_read_correct_test, single_trunc) {
  SCOPED_TRACE(m_seq.size());
  add_kmers();
  CHECK_GT(m_seq.size(), m_params.min_good_run);
  for (unsigned i = m_seq.size() - m_params.min_good_run; i < m_seq.size();
       ++i) {
    std::string err_seq = add_error(m_seq, i);
    dna_sequence expected = m_seq.substr(0, i);
    if (i < m_params.kmer_size && i >= m_seq.size() - m_params.kmer_size) {
      expected = dna_sequence();
    }
    SCOPED_TRACE(PrintToString(i) + PrintToString(err_seq) +
                 PrintToString(m_seq.size()));
    frc_output result = run_fast_read_correct(err_seq);
    EXPECT_EQ(result.corrected, expected) << " i: " << i;
    EXPECT_EQ(result.corrections, 0) << " i: " << i;
  }
}

TEST_P(fast_read_correct_test, two_errors) {
  SCOPED_TRACE(m_seq.size());
  add_kmers();
  CHECK_GT(m_seq.size(), m_params.min_good_run);

  for (unsigned i = 0; i < m_seq.size(); ++i) {
    for (unsigned j = (i + 1); j < m_seq.size(); ++j) {
      std::string err_seq = add_error(m_seq, i);
      err_seq = add_error(err_seq, j);

      dna_sequence expected = m_seq;
      unsigned expected_corrections = 2;

      if (i < m_params.kmer_size && (j - i - 1) < m_params.kmer_size &&
          j >= m_seq.size() - m_params.kmer_size) {
        // Not enough error-free sequence to start error correction with.
        expected = dna_sequence();
        expected_corrections = 0;
      } else if ((j - i - 1) < m_params.min_good_run) {
        if (i < m_params.kmer_size) {
          expected = dna_sequence();
        } else {
          expected = expected.subseq(0, i);
        }
        expected_corrections = 0;
      } else if (j >= (m_seq.size() - m_params.min_good_run)) {
        expected = expected.subseq(0, j);
        expected_corrections = 1;
      }
      frc_output result = run_fast_read_correct(err_seq);
      EXPECT_EQ(result.corrected, expected) << " i: " << i << " j: " << j
                                            << "\n";
      EXPECT_EQ(result.corrections, expected_corrections)
          << " i: " << i << " j: " << j << "\n";
    }
  }
}

TEST_P(fast_read_correct_test, three_errors) {
  add_kmers();
  CHECK_GT(m_seq.size(), m_params.min_good_run);

  parallel_for(0, m_seq.size(), [&](size_t i) {
    for (unsigned j = (i + 1); j < m_seq.size(); ++j) {
      for (unsigned k = (j + 1); k < m_seq.size(); ++k) {
        std::string err_seq = add_error(m_seq, i);
        err_seq = add_error(err_seq, j);
        err_seq = add_error(err_seq, k);

        dna_sequence expected = m_seq;
        unsigned expected_corrections = 2;

        if (i < m_params.kmer_size && (j - i - 1) < m_params.kmer_size &&
            (k - j - 1) < m_params.kmer_size) {
          // Not enough error-free sequence to start error correction with.
          expected = dna_sequence();
          expected_corrections = 0;
        } else if ((j - i - 1) < m_params.min_good_run) {
          if (i < m_params.kmer_size) {
            expected = dna_sequence();
            expected_corrections = 0;
          } else {
            expected = expected.subseq(0, i);
            expected_corrections = 0;
          }
        } else if ((k - j - 1) < m_params.min_good_run) {
          if (j < m_params.kmer_size) {
            expected = dna_sequence();
            expected_corrections = 0;
          } else {
            expected = expected.subseq(0, j);
            expected_corrections = 1;
          }
        } else {
          expected = expected.subseq(0, k);
          expected_corrections = 2;
        }
        frc_output result = run_fast_read_correct(err_seq);
        EXPECT_EQ(result.corrected, expected) << " i: " << i << " j: " << j
                                              << " k: " << k
                                              << " len: " << m_seq.size();
        EXPECT_EQ(result.corrections, expected_corrections)
            << " i: " << i << " j: " << j << " k: " << k
            << " len: " << m_seq.size();
      }
    }
  });
}

// Run the test for a bunch of different read lengths nearby multiples
// of the kmer size.
INSTANTIATE_TEST_CASE_P(
    fast_read_correct_test_N, fast_read_correct_test,
    ::testing::Combine(::testing::Values(30, 31, 32, 33, 58, 59, 60, 61, 62, 88,
                                         89, 90, 91, 92, 118, 119, 120, 121,
                                         122, 260),
                       ::testing::Values(error_behavior::REPLACE_WITH_N)));

INSTANTIATE_TEST_CASE_P(
    fast_read_correct_test_change, fast_read_correct_test,
    ::testing::Combine(::testing::Values(30, 31, 32, 33, 58, 59, 60, 61, 62, 88,
                                         89, 90, 91, 92, 118, 119, 120, 121,
                                         122, 260),
                       ::testing::Values(error_behavior::DIFFERENT_BASE)));

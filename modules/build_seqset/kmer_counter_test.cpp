#include "modules/build_seqset/kmer_counter.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/io/parallel.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace build_seqset {

using namespace testing;
using namespace dna_testutil;

struct kmer_counter_test_param {
  unsigned partitions;
  unsigned exact_passes;
  unsigned partition_batch_size;
};

class kmer_counter_test : public TestWithParam<kmer_counter_test_param> {
 public:
  using my_kmer_t = kmer_t;

  void SetUp() override {
    auto params = GetParam();
    std::cout << "Starting test with partitions: " << params.partitions
              << " exact_passes: " << params.exact_passes
              << " partition_batch_size: " << params.partition_batch_size
              << "\n";
    std::cout.flush();
    m_options.partitions = params.partitions;
    m_options.force_exact_passes = params.exact_passes;
    m_options.partition_batch_size = params.partition_batch_size;
  }

  template <typename processor_t>
  void submit_seqs() {
    std::random_device rand_dev;
    std::mt19937 rand_source(rand_dev());
    std::shuffle(m_seqs.begin(), m_seqs.end(), rand_source);
    parallel_for(0, 31, [this](size_t mod) {
      processor_t p(m_counter.get());
      for (size_t i = 0; i < m_seqs.size(); ++i) {
        if ((i % 31) != mod) {
          continue;
        }
        p.add(m_seqs[i].as_string());
      }
    });
  }

  void run() {
    m_counter.emplace(m_options);

    m_counter->start_prob_pass();
    submit_seqs<kmer_counter::prob_pass_processor>();
    m_counter->close_prob_pass();

    for (unsigned i = 0; i < m_counter->exact_passes(); ++i) {
      m_counter->start_exact_pass(i);
      submit_seqs<kmer_counter::exact_pass_processor>();
    }
    m_counter->close_exact_passes();
    m_counter->extract_exact_counts([&](kmer_counter::extract_iterator start,
                                        kmer_counter::extract_iterator limit) {
      std::lock_guard<std::mutex> l(m_mu);
      for (auto it = start; it != limit; ++it) {
        EXPECT_NE(it->kmer, std::numeric_limits<kmer_t>::max());

        EXPECT_TRUE(it->fwd_count || it->rev_count);
        EXPECT_TRUE(
            counts
                .emplace(it->kmer, std::make_pair(it->fwd_count, it->rev_count))
                .second)
            << it->kmer << ": " << it->fwd_count << ", " << it->rev_count;

        EXPECT_EQ(it->fwd_starts_read,
                  bool(m_expected_fwd_flags.count(it->kmer)));
        EXPECT_EQ(it->rev_starts_read,
                  bool(m_expected_rev_flags.count(it->kmer)));
      }
    });
    m_counter->close();
    m_counter.reset();

    EXPECT_THAT(counts, ContainerEq(expected_counts));
  }

  void save_expected() {
    for (const auto& seq : m_seqs) {
      auto it = seq.begin();
      auto end = seq.end();

      while (end - it >= m_options.kmer_size) {
        kmer_t kmer = make_kmer(it, m_options.kmer_size);
        bool flipped;
        kmer_t canon = canonicalize(kmer, m_options.kmer_size, flipped);

        bool is_start = (it == seq.begin());
        bool is_end = ((it + m_options.kmer_size) == seq.end());

        auto& e = expected_counts[canon];
        if (flipped) {
          e.second += m_multiply;
        } else {
          e.first += m_multiply;
        }

        if (flipped) {
          std::swap(is_start, is_end);
        }
        if (is_start) {
          m_expected_fwd_flags.insert(canon);
        }
        if (is_end) {
          m_expected_rev_flags.insert(canon);
        }

        ++it;
      }
    }
  }

 protected:
  std::vector<dna_sequence> m_seqs;
  count_kmer_options m_options;
  boost::optional<kmer_counter> m_counter;
  size_t m_multiply = 1;
  std::mutex m_mu;
  std::map<kmer_t, std::pair<size_t, size_t>> counts;
  std::map<kmer_t, std::pair<size_t, size_t>> expected_counts;
  std::set<kmer_t> m_expected_fwd_flags;
  std::set<kmer_t> m_expected_rev_flags;
};

TEST_P(kmer_counter_test, simple) {
  m_seqs.push_back(tseq("abcdefg"));
  m_seqs.push_back(tseq("abcdefg"));
  m_seqs.push_back(tseq("abcdefg"));
  m_seqs.push_back(tseq_rc("hij"));
  m_seqs.push_back(tseq_rc("hij"));
  m_seqs.push_back(tseq_rc("hij"));

  save_expected();
  run();
}

TEST_P(kmer_counter_test, overflow) {
  m_options.overflow_table_size_ratio = 1;

  for (size_t i = 0; i < 300; i++) {
    m_seqs.push_back(tseq("abcdefg"));
    m_seqs.push_back(tseq_rc("hij"));
  }
  save_expected();
  run();
}

TEST_P(kmer_counter_test, lots) {
  if (m_options.partitions == 1 || m_options.force_exact_passes == 1) {
    // Can't split up; make sure we can fit it all in RAM at once.
    m_options.max_memory_bytes = 100ULL * 1024 * 1024;
  }
  std::mt19937 rand_source;
  for (size_t i = 0; i < 3000; i++) {
    dna_sequence seq = rand_dna_sequence(rand_source, 200);
    m_seqs.push_back(seq);
    m_seqs.push_back(seq);
    m_seqs.push_back(seq);
  }
  save_expected();
  run();
}

TEST_P(kmer_counter_test, filtered) {
  // Only kmers that occur 3 or more times should show up in the results.
  m_seqs.push_back(tseq("three"));
  m_seqs.push_back(tseq("three"));
  m_seqs.push_back(tseq("three"));
  save_expected();
  m_seqs.push_back(tseq("TWO"));
  m_seqs.push_back(tseq("TWO"));
  m_seqs.push_back(tseq("ONE"));
  run();
}

INSTANTIATE_TEST_CASE_P(  //
    kmer_counter_test, kmer_counter_test,
    ::testing::Values(  //
        kmer_counter_test_param{
            .partitions = 1, .exact_passes = 0, .partition_batch_size = 1},
        kmer_counter_test_param{
            .partitions = 1, .exact_passes = 1, .partition_batch_size = 1},
        kmer_counter_test_param{
            .partitions = 1, .exact_passes = 1, .partition_batch_size = 100000},
        kmer_counter_test_param{
            .partitions = 64, .exact_passes = 0, .partition_batch_size = 64},
        kmer_counter_test_param{
            .partitions = 64, .exact_passes = 1, .partition_batch_size = 10000},
        kmer_counter_test_param{
            .partitions = 23, .exact_passes = 1, .partition_batch_size = 59},
        kmer_counter_test_param{
            .partitions = 23, .exact_passes = 3, .partition_batch_size = 59}));

}  // namespace build_seqset

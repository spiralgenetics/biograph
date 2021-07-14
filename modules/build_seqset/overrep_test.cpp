#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/bio_mapred/kmerize_bf.h"
#include "modules/build_seqset/kmer_counter.h"
#include "modules/io/config.h"
#include "modules/mapred/output_stream.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;
using namespace build_seqset;
using namespace dna_testutil;

using read_output_type = std::vector<std::pair<read_id, unaligned_reads>>;

enum overrep_test_type { ENABLE_OVERREP, OVERREP_THRESHOLD_TOO_HIGH, RND_ERR_THRESH_TOO_LOW };

class overrep_test : public TestWithParam<overrep_test_type /* disable overrep */> {
 public:
  void SetUp() override {
    m_test_type = GetParam();

    static size_t n = 0;

    boost::filesystem::create_directories(CONF_S(temp_root));

    m_ref_path_prefix = CONF_S(temp_root) + "/ref";
    m_ref_path_prefix += std::to_string(n);

    m_repo_path = CONF_S(temp_root) + "/repo";
    m_repo_path += std::to_string(n);

    m_count_opts.kmer_size = k_kmer_size;
    m_count_opts.max_memory_bytes = 1024 * 1024 * 1024;
    m_count_opts.min_count = k_min_kmer_count;
    m_count_opts.max_prob_table_entries = 1024 * 1024;

    m_kbf_opts.kmer_size = k_kmer_size;
    m_kbf_opts.error_rate = 0.05;
    m_kbf_opts.reference = "";
    m_kbf_opts.memory_bound = 1024 * 1024 * 1024;
    m_kbf_opts.num_threads = 1;
    m_kbf_opts.min_count = k_min_kmer_count;
    m_kbf_opts.ref_size = 1024 * 1024;

    m_kbf_opts.sys_err_thresh = 1. * k_regular_coverage / (k_overrep_coverage - 1);
    m_kbf_opts.rnd_err_thresh = 1. * k_regular_coverage / (k_overrep_coverage - 1);
    m_kbf_opts.overrep = k_overrep_coverage;
    switch (m_test_type) {
      case ENABLE_OVERREP:
        break;
      case RND_ERR_THRESH_TOO_LOW:
        m_kbf_opts.sys_err_thresh = 1. * k_regular_coverage / (k_overrep_coverage + 1);
        m_kbf_opts.rnd_err_thresh = 1. * k_regular_coverage / (k_overrep_coverage + 1);
        break;
      case OVERREP_THRESHOLD_TOO_HIGH:
        m_kbf_opts.overrep = k_overrep_coverage + 1;
        break;
    }

    ++n;
  }

  void add_expected_circular_coverage(dna_sequence seq) {
    for (const auto& read : get_circular_coverage(seq)) {
      m_reads_to_kmerize.push_back(read.as_string());

      CHECK_EQ(read.size(), k_kmer_size);
      m_expected_kmers.insert(canonicalize(read.as_kmer(), k_kmer_size));
    }
  }

  void add_circular_coverage(dna_sequence seq) {
    for (const auto& read : get_circular_coverage(seq)) {
      m_reads_to_kmerize.push_back(read.as_string());
    }
  }

  // Generates reads that are each 1 kmer long.
  std::vector<dna_sequence> get_circular_coverage(dna_sequence seq) {
    CHECK_LT(k_kmer_size, seq.size());

    std::vector<dna_sequence> result;

    dna_sequence double_seq = seq;
    double_seq += seq;

    for (unsigned pos = 0; pos < seq.size(); pos++) {
      dna_sequence read_seq = double_seq.subseq(pos, k_kmer_size);

      result.push_back(read_seq);
    }
    return result;
  }

  void kmerize() {
    output_stream_params osp;
    osp.encoding = "null";
    std::unique_ptr<kv_sink> sink = osp.build(CONF_S(path_bulkdata), "all_reads", m_reads_manifest);
    CHECK(sink);
    int n = 0;
    for (const auto& seq : m_reads_to_kmerize) {
      read_id id;
      id.pair_name = std::to_string(n);
      ++n;

      unaligned_reads reads;
      reads.emplace_back();
      auto& read = reads.back();
      read.sequence = seq;

      sink->write_msgpack(id, reads);
    }
    sink->close();
    sink.reset();

    auto kmers_and_hist = run_kmerize_subtask(m_kbf_opts, m_reads_manifest, nullptr);

    m_ks = std::move(kmers_and_hist.first);

    for (const auto& k : *m_ks) {
      EXPECT_EQ(k, canonicalize(k, k_kmer_size));
      m_actual_kmers.insert(k);
    }
  }

 protected:
  static constexpr unsigned k_underrep_coverage = 2;
  static constexpr unsigned k_regular_coverage = 3;
  static constexpr unsigned k_overrep_coverage = 40;
  static constexpr unsigned k_kmer_size = 30;
  static constexpr unsigned k_min_kmer_count = 3;

  std::string m_ref_path_prefix;
  std::string m_repo_path;

  std::vector<std::string> m_reads_to_kmerize;

  manifest m_reads_manifest;

  build_seqset::count_kmer_options m_count_opts;

  kmerize_bf_params m_kbf_opts;
  std::unique_ptr<kmer_set> m_ks;

  std::set<kmer_t> m_expected_kmers;
  std::set<kmer_t> m_actual_kmers;

  overrep_test_type m_test_type;
};

TEST_P(overrep_test, test_overrep) {
  dna_sequence regular_seq = tseq("abcdef");
  dna_sequence overrep_seq = tseq("ABCDEF");
  dna_sequence underrep_seq = tseq("012345");
  dna_sequence overrep_err_seq = overrep_seq;

  for (unsigned i = 0; i < overrep_err_seq.size(); i += k_kmer_size) {
    overrep_err_seq[i] = dna_base((int(overrep_err_seq[i]) + 1) % 4);
  }

  ASSERT_EQ(regular_seq.size() % k_kmer_size, 0);
  ASSERT_EQ(overrep_seq.size() % k_kmer_size, 0);
  ASSERT_EQ(underrep_seq.size() % k_kmer_size, 0);
  ASSERT_EQ(overrep_err_seq.size() % k_kmer_size, 0);

  for (unsigned i = 0; i < k_overrep_coverage; ++i) {
    add_expected_circular_coverage(overrep_seq);
  }
  for (unsigned i = 0; i < k_regular_coverage; ++i) {
    if (m_test_type == ENABLE_OVERREP) {
      // If overrep filtering is enabled, these should be filtered out.
      add_circular_coverage(overrep_err_seq);
    } else {
      add_expected_circular_coverage(overrep_err_seq);
    }
    add_expected_circular_coverage(regular_seq);
  }

  for (unsigned i = 0; i < k_underrep_coverage; ++i) {
    add_circular_coverage(underrep_seq);
  }

  kmerize();

  EXPECT_EQ(m_actual_kmers.size(), m_expected_kmers.size());

  std::set<dna_sequence> missing_from_actual;
  std::set<dna_sequence> missing_from_expected;
  for (const auto& k : m_actual_kmers) {
    dna_sequence s(k, k_kmer_size);
    missing_from_expected.insert(s);
  }
  for (const auto& k : m_expected_kmers) {
    dna_sequence s(k, k_kmer_size);
    if (missing_from_expected.count(s)) {
      missing_from_expected.erase(s);
    } else {
      missing_from_actual.insert(s);
    }
  }

  EXPECT_THAT(missing_from_actual, IsEmpty());
  EXPECT_THAT(missing_from_expected, IsEmpty());
}

INSTANTIATE_TEST_CASE_P(  //
    overrep_test_enabling_overrep, overrep_test, ::testing::Values(ENABLE_OVERREP));
INSTANTIATE_TEST_CASE_P(  //
    overrep_test_overrep_thresh_too_high, overrep_test,
    ::testing::Values(OVERREP_THRESHOLD_TOO_HIGH));
INSTANTIATE_TEST_CASE_P(  //
    overrep_test_rnd_err_thresh_too_low, overrep_test, ::testing::Values(RND_ERR_THRESH_TOO_LOW));

#include "modules/build_seqset/read_importer.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace testing;
using namespace build_seqset;

using output_type = std::vector<std::pair<read_id, unaligned_reads>>;

struct params {
output_type* output = nullptr;
  std::mutex* mu = nullptr;
};

class import_state : public parallel_local {
 public:
  using init_type = params;

  import_state() = delete;
  import_state(const import_state&) = delete;
  import_state(const params& val) : m_params(val) {  }
  ~import_state() { close(); }

  void process(const std::vector<std::pair<read_id, unaligned_reads>>& reads) {
    for (const auto& read : reads) {
      m_local_output.emplace_back(read);
    }
  }

 private:
  void close() {
    std::lock_guard<std::mutex> l(*m_params.mu);
    for (auto& read : m_local_output) {
      m_params.output->push_back(std::move(read));
    }
    m_local_output.clear();
  }

  params m_params;
  output_type m_local_output;
};

}  // namespace

class read_importer_test : public Test {
 public:
  read_importer_test() {
    m_params.output = &m_output;
    m_params.mu = &m_mu;
    m_importer.emplace(m_params);
  }

  void expect_has_read(const std::string& read_name, const std::string& seq) {
    size_t seen_count = 0;
    for (const auto& read : m_output) {
      if (read.first.pair_name == read_name) {
        seen_count++;
        EXPECT_EQ(1, read.second.size()) << read_name;
        if (read.second.size() >= 1) {
          EXPECT_EQ(seq, read.second[0].sequence) << read_name;
        }
      }
    }
    EXPECT_EQ(1, seen_count) << read_name;
  }

  void expect_has_paired_read(const std::string& read_name,
                              const std::string& seq,
                              const std::string& pair_seq) {
    size_t seen_count = 0;
    for (const auto& read : m_output) {
      if (read.first.pair_name == read_name) {
        seen_count++;
        EXPECT_EQ(2, read.second.size()) << read_name;
        if (read.second.size() >= 2) {
          EXPECT_EQ(seq, read.second[0].sequence) << read_name;
          EXPECT_EQ(pair_seq, read.second[1].sequence) << read_name;
        }
      }
    }
    EXPECT_EQ(1, seen_count) << read_name;
  }

  std::mutex m_mu;
  std::vector<std::pair<read_id, unaligned_reads>> m_output;
  params m_params;
  boost::optional<read_importer<import_state>> m_importer;
};

TEST_F(read_importer_test, fastq) {
  m_importer->queue_fastq("golden/E_coli_phred64.fq", "");
  m_importer->import();
  EXPECT_FALSE(m_importer->got_paired());
  // Spot check.
  expect_has_read(
      "6000:1:1101:1049:2117/1",
      "GAAACCGTTGCAGGAAACGTAACCGCGGCAGCGTCAGACACAGCCAGTTGTGTCGATTGCGGTTCCACAGGC"
      "GCTTCCACTGTGCGGCTTTTTATATATA");
  expect_has_read("6000:1:1101:1042:2228/2",
                  "CGGATGTCCGTTGGCAGTGGGTGTTTATCGGCACGGCGGTGGTCTTTTTCTTCCAGCTTT"
                  "TGCGACCGGCTTTCCAGAAAGGGTTGAAAAGCGTTTCCGG");
  expect_has_read("6000:1:1101:1436:2162/2",
                  "TTTTCAGGGCTTCTTCGCTGGCGGACGGCGCAATAATCACTTCGACAAACTGACGAGAAA"
                  "TGATGGCCTGTGCGGTTTCCGCATCCAGCTCGCGGGTAAA");
  EXPECT_EQ(10, m_output.size());
}

TEST_F(read_importer_test, pair_fastq) {
  m_importer->queue_fastq("golden/E_coli_phred64.fq", "golden/quick_e_coli.fq");
  m_importer->import();
  EXPECT_TRUE(m_importer->got_paired());
  // Spot check.
  expect_has_paired_read(
      "6000:1:1101:1049:2117/1",
      "GAAACCGTTGCAGGAAACGTAACCGCGGCAGCGTCAGACACAGCCAGTTGTGTCGATTGCGGTTCCACAGGC"
      "GCTTCCACTGTGCGGCTTTTTATATATA",
      "GGTGGCTGGTTATTCGAAGGTATGTCATGGCTGTTTATGCACCTGAACAGTAATCCGTTCGGTTGTGCGGTT"
      "TTGGCCGGGCTGTTCCTG");
  expect_has_paired_read(
      "6000:1:1101:1042:2228/2",
      "CGGATGTCCGTTGGCAGTGGGTGTTTATCGGCACGGCGGTGGTCTTTTTCTTCCAGCTTT"
      "TGCGACCGGCTTTCCAGAAAGGGTTGAAAAGCGTTTCCGG",
      "GCGCTTGTTTTTATGAAGTAAAAGAATAACGGCACTTTTTGGTGAATTTGCACTCCAAGCAACGTTATTGAA"
      "TAACCAAAGGCAGTGACA");
  expect_has_paired_read(
      "6000:1:1101:1436:2162/2",
      "TTTTCAGGGCTTCTTCGCTGGCGGACGGCGCAATAATCACTTCGACAAACTGACGAGAAA"
      "TGATGGCCTGTGCGGTTTCCGCATCCAGCTCGCGGGTAAA",
      "TGACTGGCCTCAGATTGTTGACCAAGTGCGCGTTGTACACGCCGGATGCGGCGTGAACGC"
      "CTTATCCGGCCTACGAAATCGTGCTAATTC");
  expect_has_read("r0_10",
                  "GTCCGTTTCATGATATCAGTCCAGATTGACGTTACGGCAGCCAATGAGCGTGGTGAAAGT"
                  "AAACCCGCAAACCCGTGCCACCAGAATCCC");
  expect_has_read("r0_2222",
                  "GGCAGTTTTGCGTTTGTCAGCACTCTCAGACCAGCCAGTAACATTACTGACTGGCC"
                  "TTTTTATTACTTCTGCTTTAACGCCGCATACACC");

  EXPECT_EQ(2223, m_output.size());
}

TEST_F(read_importer_test, pair_fastq_with_cut1) {
  m_importer->set_cut_region(0, 60);
  m_importer->queue_fastq("golden/E_coli_phred64.fq", "golden/quick_e_coli.fq");
  m_importer->import();
  EXPECT_TRUE(m_importer->got_paired());
  // Spot check.
  expect_has_paired_read("6000:1:1101:1049:2117/1",
                         "GAAACCGTTGCAGGAAACGTAACCGCGGCAGCGTCAGACACAGCCAGTTGTGTCGATTGC",
                         "GGTGGCTGGTTATTCGAAGGTATGTCATGGCTGTTTATGCACCTGAACAGTAATCCGTTC");
  expect_has_paired_read("6000:1:1101:1042:2228/2",
                         "CGGATGTCCGTTGGCAGTGGGTGTTTATCGGCACGGCGGTGGTCTTTTTCTTCCAGCTTT",
                         "GCGCTTGTTTTTATGAAGTAAAAGAATAACGGCACTTTTTGGTGAATTTGCACTCCAAGC");
  expect_has_read("r0_10", "GTCCGTTTCATGATATCAGTCCAGATTGACGTTACGGCAGCCAATGAGCGTGGTGAAAGT");
  expect_has_read("r0_2222", "GGCAGTTTTGCGTTTGTCAGCACTCTCAGACCAGCCAGTAACATTACTGACTGGCCTTTT");

  EXPECT_EQ(2223, m_output.size());
}

TEST_F(read_importer_test, pair_fastq_with_cut2) {
  m_importer->set_cut_region(10, 25);
  m_importer->queue_fastq("golden/E_coli_phred64.fq", "golden/quick_e_coli.fq");
  m_importer->import();
  EXPECT_TRUE(m_importer->got_paired());
  // Spot check.
  expect_has_paired_read("6000:1:1101:1049:2117/1", "CAGGAAACGTAACCG", "TATTCGAAGGTATGT");
  expect_has_paired_read("6000:1:1101:1042:2228/2", "TTGGCAGTGGGTGTT", "TTATGAAGTAAAAGA");
  expect_has_read("r0_10", "TGATATCAGTCCAGA");
  expect_has_read("r0_2222", "CGTTTGTCAGCACTC");

  EXPECT_EQ(2223, m_output.size());
}

// TODO(nils): Test importing from sources other than fastq.

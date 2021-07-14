#include "modules/bio_base/seqset.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/seqset_testutil.h"

using namespace dna_testutil;
using namespace testing;

class seqset_find_test : public Test {
 public:
  seqset_find_test() { m_seqset.emplace("golden/e_coli_merged.bg/seqset"); }

  void init_rand_seq(size_t seq_count, int prefix_step = 1) {
    size_t seed = std::random_device()();
    std::cout << "Gathering random set of entries from seqset, seed = " << seed << "...\n";
    std::mt19937_64 rand_source(seed);
    std::uniform_int_distribution<uint64_t> rand_seqset_id(0, m_seqset->size() - 1);
    for (size_t i = 0; i != seq_count; ++i) {
      uint64_t seqset_id = rand_seqset_id(rand_source);
      dna_sequence seq = m_seqset->ctx_entry(seqset_id).sequence();
      m_full_entries.insert(seq);
    }

    std::cout << "Gathering seqset prefixes...\n";
    for (const auto& full : m_full_entries) {
      dna_slice slice(full.begin(), full.end());
      for (int i = full.size(); i > 0; i -= prefix_step) {
        if (!m_prefix_entries.insert(slice.subseq(0, i)).second) {
          break;
        }
      }
    }
    std::cout << "Done\n";
  }

 protected:
  boost::optional<seqset> m_seqset;
  std::set<dna_sequence> m_full_entries;
  std::set<dna_slice> m_prefix_entries;
};

TEST_F(seqset_find_test, find_existing) {
  init_rand_seq(1000);
  for (const auto& entry : m_prefix_entries) {
    auto r = m_seqset->find(entry);
    CHECK(r.valid());
    ASSERT_TRUE(r.valid());
    EXPECT_EQ(r.begin(), m_seqset->find_existing(entry)) << "Entry: " << entry;
  }
}

TEST_F(seqset_find_test, find_existing_unique) {
  init_rand_seq(1000);
  for (const auto& entry : m_prefix_entries) {
    auto r = m_seqset->find(entry);
    ASSERT_TRUE(r.valid());
    for (size_t unique_len = 1; unique_len < entry.size(); ++unique_len) {
      EXPECT_EQ(r.begin(), m_seqset->find_existing_unique(entry, unique_len))
          << "EntrY: " << entry << " unique_len=" << unique_len;
    }
  }
}

TEST_F(seqset_find_test, shared_prefix_length) {
  init_rand_seq(20, 30);

  std::vector<std::pair<seqset_range, dna_slice>> entries;
  for (const auto& seq : m_prefix_entries) {
    entries.push_back(std::make_pair(m_seqset->find(seq), seq));
  }
  for (const auto& e1 : entries) {
    seqset_range r1 = e1.first;
    dna_slice seq1 = e1.second;
    for (const auto& e2 : entries) {
      seqset_range r2 = e2.first;
      dna_slice seq2 = e2.second;

      EXPECT_EQ(r1.shared_prefix_length(r2), seq1.shared_prefix_length(seq2));
    }
  }
}

#include "modules/build_seqset/part_repo.h"
#include "modules/build_seqset/part_counts.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace testing;
using namespace build_seqset;
using namespace dna_testutil;

}  // namespace

class part_repo_test : public TestWithParam<unsigned /* depth */> {
 public:
  void SetUp() override {
    m_depth = GetParam();

    static size_t n = 0;

    boost::filesystem::create_directories(CONF_S(temp_root));

    m_ref_path_prefix = CONF_S(temp_root) + "/ref";
    m_ref_path_prefix += std::to_string(n);

    m_repo_path = CONF_S(temp_root) + "/repo";
    m_repo_path += std::to_string(n);

    m_entries.emplace(m_depth, m_ref_path_prefix, m_repo_path);
    m_entries->open_write_pass("initial");

    ++n;
  }

  void add_seq(dna_slice seq, unsigned fwd_suffixes = 1,
               unsigned rc_suffixes = 0) {
    m_entries->write(seq, fwd_suffixes, rc_suffixes);
    auto fwd_cur = seq.begin();
    while (fwd_suffixes) {
      m_sequences.emplace_back(fwd_cur, seq.end());

      --fwd_suffixes;
      if (fwd_suffixes) {
        ASSERT_TRUE(fwd_cur != seq.end());
        ++fwd_cur;
      }
    }

    auto rev_cur = seq.rcbegin();
    while (rc_suffixes) {
      m_sequences.emplace_back(rev_cur, seq.rcend());
      --rc_suffixes;
      if (rc_suffixes) {
        ASSERT_TRUE(rev_cur != seq.rcend());
        ++rev_cur;
      }
    }
  }

  void add_expected(dna_slice seq) {
    m_sequences.emplace_back(seq.begin(), seq.end());
  }

  void load_repo() {
    m_entries->flush();
    std::sort(m_sequences.begin(), m_sequences.end());
  }

  std::vector<dna_sequence> stored_sequences() {
    CHECK(m_entries);
    std::vector<dna_sequence> seqs;
    std::mutex mu;
    m_entries->for_each_partition(  //
        "initial",                  //
        [&](const part_repo::partition_ref& part) {
          std::lock_guard<std::mutex> l(mu);
          for (const auto& e : *part.main) {
            seqs.push_back(e.sequence());
          }
        });
    std::sort(seqs.begin(), seqs.end());
    return seqs;
  }

  template <typename ref_type>
  void expect_compare(dna_slice lhs, const ref_type& rlhs, dna_slice rhs,
                      const ref_type& rrhs,
                      dna_compare_result expected_result) {
    ASSERT_EQ(expected_result, lhs.compare_to(rhs)) << "lhs: " << lhs
                                                    << " rhs: " << rhs;
    EXPECT_EQ(expected_result, rlhs.compare_to(rrhs)) << "lhs: " << lhs
                                                      << " rhs: " << rhs;
    EXPECT_EQ(expected_result, rlhs.sequence().compare_to(rhs))
        << "lhs: " << lhs << " rhs: " << rhs;
    unsigned expected_shared = 0;
    while (expected_shared < lhs.size() && expected_shared < rhs.size() &&
           lhs[expected_shared] == rhs[expected_shared]) {
      expected_shared++;
    }

    bool expect_lt = false;
    bool converse_expect_lt = false;
    dna_compare_result converse_expected = dna_compare_result::EQUAL;
    switch (expected_result) {
      case dna_compare_result::FIRST_IS_LESS:
        converse_expected = dna_compare_result::SECOND_IS_LESS;
        expect_lt = true;
        converse_expect_lt = false;
        break;
      case dna_compare_result::FIRST_IS_PREFIX:
        converse_expected = dna_compare_result::SECOND_IS_PREFIX;
        expect_lt = true;
        converse_expect_lt = false;
        break;
      case dna_compare_result::EQUAL:
        converse_expected = dna_compare_result::EQUAL;
        expect_lt = false;
        converse_expect_lt = false;
        break;
      case dna_compare_result::SECOND_IS_PREFIX:
        converse_expected = dna_compare_result::FIRST_IS_PREFIX;
        expect_lt = false;
        converse_expect_lt = true;
        break;
      case dna_compare_result::SECOND_IS_LESS:
        converse_expected = dna_compare_result::SECOND_IS_LESS;
        expect_lt = false;
        converse_expect_lt = true;
        break;
    }
    ASSERT_EQ(converse_expected, rhs.compare_to(lhs));
    EXPECT_EQ(converse_expected, rrhs.compare_to(rlhs));
    EXPECT_EQ(converse_expected, rrhs.sequence().compare_to(lhs));

    EXPECT_EQ(expect_lt, rlhs < rrhs);
    EXPECT_EQ(expect_lt, rrhs > rlhs);

    EXPECT_EQ(converse_expect_lt, rrhs < rlhs);
    EXPECT_EQ(converse_expect_lt, rlhs > rrhs);

    ASSERT_EQ(expected_shared, lhs.shared_prefix_length(rhs));
    EXPECT_EQ(expected_shared, rlhs.shared_prefix_length(rrhs));
  }

 protected:
  std::string m_ref_path_prefix;
  std::string m_repo_path;

  boost::optional<part_repo> m_entries;
  std::vector<dna_sequence> m_sequences;

  unsigned m_depth;
};

TEST_P(part_repo_test, fwd_and_rev_count_simple) {
  add_seq(tseq("abcde"), 0, 1);
  add_seq(tseq("fghij"), 1, 1);
  load_repo();

  EXPECT_THAT(stored_sequences(), ElementsAreArray(m_sequences));
}

TEST_P(part_repo_test, fwd_and_rev_count) {
  add_seq(tseq("abcde"), 0, 1);
  add_seq(tseq("fghij"), tseq("fghij").size(), 0);
  add_seq(tseq("klmno"), 0, tseq("klmno").size());
  add_seq(tseq("pqrstu"), tseq("pqrstu").size() / 3, tseq("pqrstu").size() / 2);
  load_repo();

  EXPECT_THAT(stored_sequences(), ElementsAreArray(m_sequences));
}

TEST_P(part_repo_test, sort) {
  add_seq(tseq("abcde"), 0, 1);
  add_seq(tseq("fghij"), tseq("fghij").size(), 0);
  add_seq(tseq("klmno"), 0, tseq("klmno").size());
  add_seq(tseq("pqrstu"), tseq("pqrstu").size() / 3, tseq("pqrstu").size() / 2);

  m_entries->flush();
  m_entries->for_each_partition(
      "initial", [&](const part_repo::partition_ref& part) {
        std::vector<seq_repository::entry_data> data(part.main->data_begin(),
                                                     part.main->data_end());
        std::sort(data.begin(), data.end(),
                  part.main->less_than_using_repo());
        auto b = m_entries->open_ref_builder(part.part_id, "sorted");
        for (const auto& e : data) {
          b->write_entry(e);
        }
      });
  m_entries->flush();
  std::mutex mu;
  m_entries->for_each_partition(
      "sorted", [&](const part_repo::partition_ref& part) {
        std::lock_guard<std::mutex> l(mu);
        auto it = part.main->begin();
        if (it == part.main->end()) {
          return;
        }
        auto next = it;
        ++next;
        while (next != part.main->end()) {
          EXPECT_LE(*it, *next) << it->sequence() << " <= " << next->sequence();
          it = next;
          ++next;
        }
      });
}

TEST_P(part_repo_test, pushes) {
  add_seq(tseq("abcde"), 0, 1);
  // add_seq(tseq("fghij"), tseq("fghij").size(), 0);
  // add_seq(tseq("klmno"), 0, tseq("klmno").size());
  // add_seq(tseq("pqrstu"), tseq("pqrstu").size() / 3, tseq("pqrstu").size() /
  // 2);

  m_entries->flush();

  std::map<dna_sequence, dna_base_array<std::set<dna_sequence>>>
      expected_pushes;
  for (const auto& seq : m_sequences) {
    if (!seq.size()) {
      continue;
    }
    dna_sequence popped_seq = seq.subseq(1, seq.size() - 1);
    dna_sequence pref = popped_seq;
    while (pref.size() < m_depth) {
      pref.push_back(dna_base('A'));
    }
    expected_pushes[pref.subseq(0, m_depth)][seq[0]].insert(seq);
  }

  std::mutex mu;
  dna_base_array<std::set<dna_sequence>> actual_pushes;
  m_entries->for_each_partition(
      "initial", [&](const part_repo::partition_ref& part) {
        std::lock_guard<std::mutex> l(mu);

        std::cout << "Partition " << part.prefix << " with "
                  << part.main->end() - part.main->begin() << " entries:\n";
        for (dna_base b : dna_bases()) {
          std::set<dna_sequence> actual;
          for (auto it = part.pushed[b].first; it != part.pushed[b].second;
               ++it) {
            actual.insert(it->sequence());
          }
          EXPECT_EQ(actual, expected_pushes[part.prefix][b])
              << "Prefix " << part.prefix << " base " << b;
        }
      });
}

INSTANTIATE_TEST_CASE_P(  //
    part_repo_test, part_repo_test, ::testing::Values(1, 2, 3, 4));

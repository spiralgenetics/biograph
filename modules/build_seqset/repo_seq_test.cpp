#include "modules/build_seqset/repo_seq.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/io/config.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace testing;
using namespace build_seqset;
using namespace dna_testutil;

}  // namespace

TEST(repo_seq_self_contained_test, self_contained) {
  dna_sequence seq;

  for (unsigned len = 0; len < seq_repository::k_inline_bases * 2; ++len) {
    unsigned inline_len =
        std::min<unsigned>(len, seq_repository::k_inline_bases);
    seq_repository::entry_data edata(len, seq.subseq(0, inline_len),
                                     seq_repository::k_inline_bases, false);
    seq_repository::entry e(edata, seq, 0 /* popped */);

    EXPECT_EQ(e.sequence(), seq);

    if (len > 0) {
      seq_repository::entry popped = e.pop_front();

      EXPECT_EQ(popped.sequence(), seq.subseq(1, seq.size() - 1));
    }

    seq.push_back((__builtin_popcountl(len) & 1) ? dna_base('C')
                                                 : dna_base('T'));
  }
}

class repo_seq_test : public Test {
 public:
  void SetUp() override {
    static size_t n = 0;

    boost::filesystem::create_directories(CONF_S(temp_root));

    m_ref_path = CONF_S(temp_root) + "/ref";
    m_ref_path += std::to_string(n);
    m_ref_builder.emplace(m_ref_path);

    m_repo_path = CONF_S(temp_root) + "/repo";
    m_repo_path += std::to_string(n);
    m_repo_builder.emplace(m_repo_path);

    ++n;
  }

  void add_seq(dna_slice seq, unsigned fwd_suffixes = 1,
               unsigned rc_suffixes = 0) {
    m_ref_builder->write_sequence(seq, *m_repo_builder, fwd_suffixes,
                                  rc_suffixes);
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
    m_ref_builder.reset();
    m_repo_builder.reset();

    m_entries.reset();
    m_entries.emplace(m_ref_path, m_repo_path);
  }

  std::vector<dna_sequence> stored_sequences() {
    CHECK(m_entries);
    std::vector<dna_sequence> seqs;
    for (const auto& e : *m_entries) {
      seqs.push_back(e.sequence());
    }
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
  std::string m_ref_path;
  boost::optional<seq_repository::ref_builder> m_ref_builder;
  std::string m_repo_path;
  boost::optional<seq_repository::repo_builder> m_repo_builder;

  boost::optional<seq_repository> m_entries;
  std::vector<dna_sequence> m_sequences;
};

TEST_F(repo_seq_test, compares) {
  dna_sequence seq;
  for (unsigned len = 0; len < seq_repository::k_inline_bases * 2; ++len) {
    add_seq(seq);

    dna_sequence with_a = seq;
    with_a.push_back(dna_base('A'));
    add_seq(with_a);
    with_a.push_back(dna_base('A'));
    add_seq(with_a);

    dna_sequence with_t = seq;
    with_t.push_back(dna_base('T'));
    add_seq(with_t);
    with_t.push_back(dna_base('T'));
    add_seq(with_t);

    seq.push_back(__builtin_popcountl(len) ? dna_base('C') : dna_base('G'));
  }

  load_repo();

  auto it = m_entries->begin();
  auto sit = m_sequences.begin();
  while (it != m_entries->end()) {
    ASSERT_TRUE(it != m_entries->end() && sit != m_sequences.end());
    seq_repository::reference rseq = *it++;
    dna_sequence seq = *sit++;
    ASSERT_TRUE(it != m_entries->end() && sit != m_sequences.end());
    seq_repository::reference rseq_A = *it++;
    dna_sequence seq_A = *sit++;
    ASSERT_TRUE(it != m_entries->end() && sit != m_sequences.end());
    seq_repository::reference rseq_AA = *it++;
    dna_sequence seq_AA = *sit++;
    ASSERT_TRUE(it != m_entries->end() && sit != m_sequences.end());
    seq_repository::reference rseq_T = *it++;
    dna_sequence seq_T = *sit++;
    ASSERT_TRUE(it != m_entries->end() && sit != m_sequences.end());
    seq_repository::reference rseq_TT = *it++;
    dna_sequence seq_TT = *sit++;

    expect_compare(seq, rseq, seq, rseq, dna_compare_result::EQUAL);

    expect_compare(seq, rseq, seq_A, rseq_A,
                   dna_compare_result::FIRST_IS_PREFIX);
    expect_compare(seq, rseq, seq_AA, rseq_AA,
                   dna_compare_result::FIRST_IS_PREFIX);
    expect_compare(seq, rseq, seq_T, rseq_T,
                   dna_compare_result::FIRST_IS_PREFIX);
    expect_compare(seq, rseq, seq_TT, rseq_TT,
                   dna_compare_result::FIRST_IS_PREFIX);

    expect_compare(seq_A, rseq_A, seq_T, rseq_T,
                   dna_compare_result::FIRST_IS_LESS);
    expect_compare(seq_AA, rseq_AA, seq_T, rseq_T,
                   dna_compare_result::FIRST_IS_LESS);
    expect_compare(seq_A, rseq_A, seq_TT, rseq_TT,
                   dna_compare_result::FIRST_IS_LESS);
    expect_compare(seq_AA, rseq_AA, seq_TT, rseq_TT,
                   dna_compare_result::FIRST_IS_LESS);
  }
}

TEST_F(repo_seq_test, pop_front) {
  dna_sequence seq = tseq("abcdefg");
  add_seq(seq);
  load_repo();

  while (seq.size() != 0) {
    {
      seq_repository::ref_builder appender(m_ref_path);

      auto it = m_entries->begin();
      ASSERT_TRUE(it != m_entries->end());
      seq_repository::entry e;
      do {
        e = *it;
        ++it;
      } while (it != m_entries->end());

      EXPECT_EQ(seq, e.sequence());

      seq = seq.subseq(1, seq.size() - 1);
      add_expected(seq);
      appender.write_entry(e.pop_front());
    }
    load_repo();
  }

  auto it = m_entries->begin();
  ASSERT_TRUE(it != m_entries->end());
  seq_repository::entry e;
  do {
    e = *it;
    ++it;
  } while (it != m_entries->end());

  EXPECT_EQ(seq, e.sequence());
}

TEST_F(repo_seq_test, sorts_appends_pops) {
  dna_sequence seq;
  for (unsigned len = 0; len < seq_repository::k_inline_bases * 2; ++len) {
    add_seq(seq);

    dna_sequence with_a = seq;
    with_a.push_back(dna_base('A'));
    add_seq(with_a);
    with_a.push_back(dna_base('A'));
    add_seq(with_a);

    dna_sequence with_t = seq;
    with_t.push_back(dna_base('T'));
    add_seq(with_t);
    with_t.push_back(dna_base('T'));
    add_seq(with_t);

    seq.push_back((__builtin_popcountl(len) & 1) ? dna_base('C')
                                                 : dna_base('G'));
  }

  load_repo();

  EXPECT_EQ(m_entries->data_end() - m_entries->data_begin(),
            m_entries->end() - m_entries->begin());

  // Sort, make sure elements are in order.
  std::vector<seq_repository::entry_data> data(m_entries->data_begin(),
                                               m_entries->data_end());
  std::sort(data.begin(), data.end(), m_entries->less_than_using_repo());

  auto it = seq_repository::iterator(data.data(), m_entries->repo());
  auto end = it + data.size();
  ASSERT_TRUE(it != end);
  auto next = it;
  ++next;
  while (next != end) {
    EXPECT_LE(*it, *next) << it->sequence() << " <= " << next->sequence();
    EXPECT_LE(it->sequence(), next->sequence());

    it = next;
    ++next;
  }

  size_t entry_count = data.size();

  // Add pop front of each entry.
  size_t popped_entry_count = 0;
  {
    seq_repository::ref_builder appender(m_ref_path);

    for (const auto& e_data : data) {
      seq_repository::reference e(&e_data, m_entries->repo());
      if (e.size() > 0) {
        dna_sequence e_seq = e.sequence();
        if (popped_entry_count & 1) {
          seq_repository::entry popped = e.pop_front();
          appender.write_entry(popped);
          dna_sequence popped_seq = popped.sequence();
          EXPECT_EQ(popped_seq, e_seq.subseq(1, e_seq.size() - 1)) << "e_seq: "
                                                                   << e_seq;
          add_expected(popped_seq);
        } else {
          seq_repository::entry es = e;
          seq_repository::entry popped = es.pop_front();
          appender.write_entry(popped);
          dna_sequence popped_seq = popped.sequence();
          EXPECT_EQ(popped_seq, e_seq.subseq(1, e_seq.size() - 1)) << "e_seq: "
                                                                   << e_seq;
          add_expected(popped_seq);
        }
        ++popped_entry_count;
      }
    }
  }

  load_repo();

  data.clear();
  std::copy(m_entries->data_begin(), m_entries->data_end(),
            std::back_inserter(data));
  std::sort(data.begin(), data.end(), m_entries->less_than_using_repo());

  it = seq_repository::iterator(data.data(), m_entries->repo());
  end = it + data.size();
  ASSERT_TRUE(it != end);
  next = it;
  ++next;
  while (next != end) {
    EXPECT_LE(*it, *next) << it->sequence() << " <= " << next->sequence();
    EXPECT_LE(it->sequence(), next->sequence());

    it = next;
    ++next;
  }

  EXPECT_EQ(entry_count + popped_entry_count,
            m_entries->data_end() - m_entries->data_begin())
      << " Entries: " << entry_count
      << " Popped entries added: " << popped_entry_count;

  std::sort(m_sequences.begin(), m_sequences.end());
  auto stored = stored_sequences();
  std::sort(stored.begin(), stored.end());
  EXPECT_EQ(stored.size(), entry_count + popped_entry_count);
  EXPECT_EQ(m_sequences.size(), entry_count + popped_entry_count);
  EXPECT_THAT(stored, ElementsAreArray(m_sequences));
}

TEST_F(repo_seq_test, fwd_and_rev_count_simple) {
  add_seq(tseq("abcde"), 0, 1);
  add_seq(tseq("fghij"), 1, 1);
  load_repo();

  EXPECT_THAT(stored_sequences(), ElementsAreArray(m_sequences));
}

TEST_F(repo_seq_test, fwd_and_rev_count) {
  add_seq(tseq("abcde"), 0, 1);
  add_seq(tseq("fghij"), tseq("fghij").size(), 0);
  add_seq(tseq("klmno"), 0, tseq("klmno").size());
  add_seq(tseq("pqrstu"), tseq("pqrstu").size() / 3, tseq("pqrstu").size() / 2);
  load_repo();

  EXPECT_THAT(stored_sequences(), ElementsAreArray(m_sequences));
}

TEST_F(repo_seq_test, rc_pop_front) {
  dna_sequence seq = tseq("abcdefg");
  add_seq(seq.rev_comp(), 0, 1);
  load_repo();

  while (seq.size() != 0) {
    {
      seq_repository::ref_builder appender(m_ref_path);

      auto it = m_entries->begin();
      ASSERT_TRUE(it != m_entries->end());
      seq_repository::entry e;
      do {
        e = *it;
        ++it;
      } while (it != m_entries->end());

      EXPECT_EQ(seq, e.sequence());

      seq = seq.subseq(1, seq.size() - 1);
      add_expected(seq);
      appender.write_entry(e.pop_front());
    }
    load_repo();
  }

  auto it = m_entries->begin();
  ASSERT_TRUE(it != m_entries->end());
  seq_repository::entry e;
  do {
    e = *it;
    ++it;
  } while (it != m_entries->end());

  EXPECT_EQ(seq, e.sequence());
}

TEST_F(repo_seq_test, pop_iterator) {
  dna_sequence seq = tseq("abcde");
  add_seq(seq);

  load_repo();

  {
    seq_repository::ref_builder appender(m_ref_path);
    seq_repository::popped_iterator begin = m_entries->begin().pop_front();
    seq_repository::popped_iterator end = m_entries->end().pop_front();
    for (unsigned pops = 1; pops < 100; pops++) {
      seq_repository::entry e;
      for (auto it = begin; it != end; ++it) {
        e = *it;
        appender.write_entry(e);
      }
      if (seq.size() > 0) {
        seq = seq.subseq(1, seq.size() - 1);
      }
      add_expected(seq);
      EXPECT_EQ(e.sequence(), seq);

      begin = begin.pop_front();
      end = end.pop_front();
    }
  }

  load_repo();

  EXPECT_THAT(stored_sequences(), ElementsAreArray(m_sequences));
}

TEST_F(repo_seq_test, empty) {
  load_repo();

  EXPECT_TRUE(m_entries->begin() == m_entries->end());
  EXPECT_TRUE(m_entries->begin().pop_front() == m_entries->end().pop_front());
}

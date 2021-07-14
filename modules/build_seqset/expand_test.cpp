#include "modules/build_seqset/expand.h"
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

class expand_test : public TestWithParam<unsigned /* partition depth */> {
 public:
  void SetUp() override {
    static size_t n = 0;

    m_depth = GetParam();

    boost::filesystem::create_directories(CONF_S(temp_root));

    m_ref_path = CONF_S(temp_root) + "/ref";
    m_ref_path += std::to_string(n);
    m_entries_path = CONF_S(temp_root) + "/repo";
    m_entries_path += std::to_string(n);

    m_entries.emplace(m_depth, m_ref_path, m_entries_path);
    m_entries->open_write_pass("initial");

    ++n;
  }

  void add_seq(dna_slice seq, unsigned fwd_suffixes = 1,
               unsigned rc_suffixes = 0) {
    m_entries->write(seq, fwd_suffixes, rc_suffixes);

    for (auto it = seq.begin(); it != seq.end(); ++it) {
      add_expected_seq(dna_sequence(it, seq.end()));
    }
    if (rc_suffixes) {
      for (auto it = seq.rcbegin(); it != seq.rcend(); ++it) {
        add_expected_seq(dna_sequence(it, seq.rcend()));
      }
    }
  }

  void add_expected_seq(const dna_sequence& seq) {
    m_expected_seqs.insert(seq);
  }

  std::set<dna_sequence> dedup(std::set<dna_sequence> in) {
    std::set<dna_sequence> out;

    for (auto it = in.begin(); it != in.end(); ++it) {
      if (out.begin() != out.end()) {
        auto last_out = out.end();
        --last_out;
        dna_compare_result cmp = last_out->compare_to(*it);
        if (cmp == dna_compare_result::EQUAL ||
            cmp == dna_compare_result::FIRST_IS_PREFIX) {
          out.erase(last_out);
        } else {
          CHECK_EQ(cmp, dna_compare_result::FIRST_IS_LESS);
        }
      }
      out.insert(*it);
    }
    return out;
  }

  void run_expand() {
    m_entries->flush();

    m_expander.emplace(*m_entries, true /* save temporary files */);
    m_expander->sort_and_dedup("", "initial", "init_sorted", "", 0, 0);
    m_expander->expand("init_sorted", "init_expanded", 16, 255);
    m_expander->sort_and_dedup("init_sorted", "init_expanded", "pass2_sorted",
                               "pass2_expanded", 1, 15);
    m_expander->sort_and_dedup("pass2_sorted", "pass2_expanded", "complete", "",
                               0, 0);
    size_t more_expand_needed =
        m_expander->expand("complete", "complete_expanded", 1, 255);
    size_t dedupped = m_expander->sort_and_dedup(
        "complete", "complete_expanded", "test_out", "", 0, 0);
    EXPECT_EQ(more_expand_needed, dedupped);
    m_expander.reset();
  }

  std::vector<dna_sequence> stored_entries() {
    return pass_entries("complete");
  }

  void dump_pass_entries(const std::string& pass_name) {
    std::cout << "Pass " << pass_name << ":\n";
    for (const auto& seq : pass_entries(pass_name)) {
      std::cout << "  " << seq << "\n";
    }
    std::cout.flush();
  }

  std::vector<dna_sequence> pass_entries(const std::string& pass_name) {
    CHECK(m_entries);
    std::map<dna_sequence, std::vector<dna_sequence>> part_seqs;
    std::mutex mu;
    m_entries->for_each_partition(
        pass_name, [&](const part_repo::partition_ref& part) {
          std::lock_guard<std::mutex> l(mu);
          for (const auto& e : *part.main) {
            part_seqs[part.prefix].push_back(e.sequence());
          }
        });
    std::vector<dna_sequence> seqs;
    for (const auto& part : part_seqs) {
      std::copy(part.second.begin(), part.second.end(),
                std::back_inserter(seqs));
    }
    return seqs;
  }

  unsigned m_depth;
  std::string m_ref_path;
  std::string m_entries_path;
  boost::optional<part_repo> m_entries;

  std::vector<dna_sequence> m_sequences;
  boost::optional<expander> m_expander;

  std::set<dna_sequence> m_expected_seqs;
};

TEST_P(expand_test, sort_and_dedup) {
  add_seq(tseq("bcd"));
  add_seq(tseq("ab"));
  add_seq(tseq("a"));
  add_seq(tseq("abc"));
  add_seq(dna_sequence("C"));

  m_entries->flush();
  m_entries->open_write_pass("initial2");

  add_seq(tseq("bc"));
  add_seq(tseq("bcd"));
  add_seq(tseq("abcd"));

  m_entries->flush();
  m_expander.emplace(*m_entries, true /* keep temporary files */);
  m_expander->sort_and_dedup("", "initial", "sorted", "expanded", 0, 0);
  EXPECT_THAT(pass_entries("sorted"),
              UnorderedElementsAre(tseq("abc"), tseq("bcd")));
  EXPECT_THAT(pass_entries("expanded"), IsEmpty());

  m_expander->sort_and_dedup("sorted", "initial2", "sorted2", "expanded2", 2,
                             3);
  EXPECT_THAT(pass_entries("sorted2"),
              UnorderedElementsAre(tseq("abcd"), tseq("bcd")));
  EXPECT_THAT(pass_entries("expanded2"),
              UnorderedElementsAre(drop_front(1, tseq("abcd")),
                                   drop_front(3, tseq("abcd")),
                                   drop_front(5, tseq("abcd"))));
}

TEST_P(expand_test, sort_expand) {
  add_seq(tseq("abc"));

  m_entries->flush();
  m_entries->open_write_pass("initial2");

  add_seq(dna_sequence("T") + tseq("bcdefg"));

  m_entries->flush();

  m_expander.emplace(*m_entries, true /* keep temporary files */);
  m_expander->sort_and_dedup("initial", "initial2", "sorted", "expanded",
                             k_dna_test_sequence_length, 2);
  EXPECT_THAT(
      pass_entries("sorted"),
      UnorderedElementsAre(tseq("abc"), dna_sequence("T") + tseq("bcdefg")));
  EXPECT_THAT(pass_entries("expanded"),
              UnorderedElementsAre(tseq("bcdefg"), tseq("cdefg")));
}

TEST_P(expand_test, simple_small) {
  add_seq(tseq("a"));
  add_seq(tseq_rc("a"));

  run_expand();

  EXPECT_THAT(stored_entries(), ElementsAreArray(dedup(m_expected_seqs)));
  EXPECT_THAT(
      stored_entries(),
      ElementsAreArray({dna_sequence("AAAATTAC"), dna_sequence("AAATTAC"),
                        dna_sequence("AATTAC"), dna_sequence("AATTTTAG"),
                        dna_sequence("AC"), dna_sequence("AG"),
                        dna_sequence("ATTAC"), dna_sequence("ATTTTAG"),
                        dna_sequence("CTAAAATTAC"), dna_sequence("GTAATTTTAG"),
                        dna_sequence("TAAAATTAC"), dna_sequence("TAATTTTAG"),
                        dna_sequence("TAC"), dna_sequence("TAG"),
                        dna_sequence("TTAC"), dna_sequence("TTAG"),
                        dna_sequence("TTTAG"), dna_sequence("TTTTAG")}));
}

TEST_P(expand_test, simple_big) {
  add_seq(tseq("abcdefg"));

  run_expand();

  EXPECT_THAT(stored_entries(), ElementsAreArray(dedup(m_expected_seqs)));
}

TEST_P(expand_test, preexpanded) {
  auto seq = dna_sequence("GATTACA");
  add_seq(seq, seq.size(), seq.size());

  run_expand();

  EXPECT_THAT(stored_entries(), ElementsAreArray(dedup(m_expected_seqs)));
}

TEST_P(expand_test, bigger_preexpanded) {
  auto seq = tseq("abcdefg");
  add_seq(seq, seq.size(), seq.size());

  run_expand();

  EXPECT_THAT(stored_entries(), ElementsAreArray(dedup(m_expected_seqs)));
}

INSTANTIATE_TEST_CASE_P(  //
    expand_test, expand_test, ::testing::Values(1, 2, 3, 4));

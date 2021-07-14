#include "modules/build_seqset/correct_reads.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/build_seqset/part_repo.h"
#include "modules/io/config.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace testing;
using namespace build_seqset;
using namespace dna_testutil;

}  // namespace

class correct_reads_test : public Test {
 public:
  void SetUp() override {
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

  void set_initial_repo(dna_slice seq) { m_entries->add_initial_repo(seq); }

  void add_kmers(dna_slice seq) {
    CHECK_GE(seq.size(), m_kmer_size);
    for (auto it = seq.begin(); (it + m_kmer_size - 1) != seq.end(); ++it) {
      m_kmers.insert(canonicalize(make_kmer(it, m_kmer_size), m_kmer_size));
    }
    m_first_kmers.insert(make_kmer(seq.begin(), m_kmer_size));
    m_first_kmers.insert(make_kmer(seq.rcbegin(), m_kmer_size));
    m_kmer_bases += seq.size();
  }

  void start_correction() {
    m_ks.emplace(m_kmer_bases, m_kmer_size, get_maximum_mem_bytes(),
                 [&](const kmer_set::kmer_output_f& out_f, progress_handler_t) {
                   for (kmer_t k : m_kmers) {
                     unsigned flags = 0;
                     if (m_first_kmers.count(k)) {
                       flags |= kmer_set::k_fwd_starts_read;
                     }
                     if (m_first_kmers.count(rev_comp(k, m_kmer_size))) {
                       flags |= kmer_set::k_rev_starts_read;
                     }

                     out_f(k, flags);
                   }
                 });
    m_entries->flush();
    m_entries->open_write_pass("initial");
    m_cr.emplace(*m_entries, *m_ks, m_params);
    m_cr->add_initial_repo();
  }

  void add_read(dna_slice seq) {
    unaligned_read r;
    r.sequence = seq.as_string();
    corrected_read c;
    m_cr->correct(r, c);
  }

  void load_repo() {
    m_cr.reset();
    m_entries->flush();
  }

  std::set<dna_sequence> stored_sequences() {
    CHECK(m_entries);
    std::set<dna_sequence> seqs;
    std::mutex mu;
    m_entries->for_each_partition("initial",
                                  [&](const part_repo::partition_ref& part) {
                                    std::lock_guard<std::mutex> l(mu);
                                    for (const auto& e : *part.main) {
                                      seqs.insert(e.sequence());
                                    }
                                  });
    return seqs;
  }

  Matcher<const std::set<dna_sequence>&> ContainsAllOf(const std::set<dna_sequence>& seqs) {
    CHECK_GT(seqs.size(), 0);
    Matcher<const std::set<dna_sequence>&> result = _;
    for (const auto& seq : seqs) {
      result = AllOf(result, Contains(seq));
    }
    return result;
  }

 protected:
  std::string m_ref_path_prefix;
  std::string m_repo_path;

  size_t m_kmer_bases = 0;
  std::set<kmer_t> m_kmers;
  std::set<kmer_t> m_first_kmers;
  boost::optional<kmer_set> m_ks;
  boost::optional<part_repo> m_entries;
  boost::optional<correct_reads> m_cr;

  read_correction_params m_params;

  unsigned m_depth = 1;
  unsigned m_kmer_size = 2 * k_dna_test_sequence_length;
};

TEST_F(correct_reads_test, simple_fwd) {
  add_kmers(tseq("bcde"));
  add_kmers(dna_sequence("GA") + tseq("bcde") + dna_sequence("A"));
  start_correction();
  add_read(dna_sequence("GA") + tseq("bcde") + dna_sequence("A"));
  load_repo();
  std::set<dna_sequence> expected = {
      dna_sequence("GA") + tseq("bcde") + dna_sequence("A"),
      dna_sequence("A") + tseq("bcde") + dna_sequence("A"),
      (dna_sequence("GA") + tseq("bcde") + dna_sequence("A")).rev_comp()};

  EXPECT_THAT(stored_sequences(), ContainerEq(expected));
  EXPECT_EQ(m_entries->repo_slice().size() / 4,
            ((dna_sequence("GA") + tseq("bcde") + dna_sequence("A")).size() + 3) / 4);
}

TEST_F(correct_reads_test, simple_rc) {
  add_kmers(tseq("bcde"));
  add_kmers(dna_sequence("A") + tseq("bcde") + dna_sequence("GA"));
  start_correction();
  add_read(dna_sequence("A") + tseq("bcde") + dna_sequence("GA"));
  load_repo();
  std::set<dna_sequence> expected = {
      (dna_sequence("A") + tseq("bcde") + dna_sequence("GA")).rev_comp(),
      (dna_sequence("A") + tseq("bcde") + dna_sequence("G")).rev_comp(),
      (dna_sequence("A") + tseq("bcde") + dna_sequence("GA"))};

  EXPECT_THAT(stored_sequences(), ContainerEq(expected));
  EXPECT_EQ(m_entries->repo_slice().size() / 4,
            ((dna_sequence("A") + tseq("bcde") + dna_sequence("GA")).size() + 3) / 4);
}

TEST_F(correct_reads_test, fwd_in_repo) {
  set_initial_repo(tseq("abcdefghij"));
  add_kmers(tseq("abcd"));
  add_kmers(tseq("ghij"));
  start_correction();
  add_read(tseq("abcd"));
  add_read(tseq("ghij"));
  load_repo();
  std::set<dna_sequence> expected = {dna_sequence(tseq("abcd")), dna_sequence(tseq_rc("abcd")),
                                     dna_sequence(tseq("ghij")), dna_sequence(tseq_rc("ghij"))};

  EXPECT_THAT(stored_sequences(), ContainsAllOf(expected));
  EXPECT_EQ(m_entries->repo_slice().size() / 4, (tseq("abcdefghij").size() + 3) / 4);
}

TEST_F(correct_reads_test, rc_in_repo) {
  set_initial_repo(tseq_rc("abcdefghij"));
  add_kmers(tseq("abcd"));
  add_kmers(tseq("ghij"));
  start_correction();
  add_read(tseq("abcd"));
  add_read(tseq("ghij"));
  load_repo();
  std::set<dna_sequence> expected = {
      dna_sequence(tseq("abcd")), dna_sequence(tseq_rc("abcd")),
      dna_sequence(tseq("ghij")), dna_sequence(tseq_rc("ghij"))};

  EXPECT_THAT(stored_sequences(), ContainsAllOf(expected));
  EXPECT_EQ(m_entries->repo_slice().size() / 4,
            (tseq_rc("abcdefghij").size() + 3) / 4);
}

TEST_F(correct_reads_test, fwd_almost_in_repo) {
  set_initial_repo(tseq("abcdefghij"));
  add_kmers(dna_G + tseq("abcd"));
  add_kmers((dna_G + tseq("abcd")).subseq(1, tseq("abcd").size() - 1));
  add_kmers(tseq("ghij") + dna_G);
  add_kmers((tseq("ghij") + dna_G).subseq(1, tseq("ghij").size() - 1));
  start_correction();
  add_read(dna_G + tseq("abcd"));
  add_read(tseq("ghij") + dna_G);
  load_repo();
  std::set<dna_sequence> expected = {dna_sequence(dna_G + tseq("abcd")),
                                     dna_sequence(tseq_rc("abcd") + dna_C),
                                     dna_sequence(tseq("ghij") + dna_G),
                                     dna_sequence(dna_C + tseq_rc("ghij"))};

  EXPECT_THAT(stored_sequences(), ContainerEq(expected));
  EXPECT_GT(m_entries->repo_slice().size() / 4,
            (tseq("abcdefghij").size() + 3) / 4);
}

TEST_F(correct_reads_test, rev_almost_in_repo) {
  set_initial_repo(tseq_rc("abcdefghij"));
  add_kmers(dna_G + tseq("abcd"));
  add_kmers((dna_G + tseq("abcd")).subseq(1, tseq("abcd").size() - 1));
  add_kmers(tseq("ghij") + dna_G);
  add_kmers((tseq("ghij") + dna_G).subseq(1, tseq("ghij").size() - 1));
  start_correction();
  add_read(dna_G + tseq("abcd"));
  add_read(tseq("ghij") + dna_G);
  load_repo();
  std::set<dna_sequence> expected = {dna_sequence(dna_G + tseq("abcd")),
                                     dna_sequence(tseq_rc("abcd") + dna_C),
                                     dna_sequence(tseq("ghij") + dna_G),
                                     dna_sequence(dna_C + tseq_rc("ghij"))};

  EXPECT_THAT(stored_sequences(), ContainerEq(expected));
  EXPECT_GT(m_entries->repo_slice().size() / 4,
            (tseq("abcdefghij").size() + 3) / 4);
}

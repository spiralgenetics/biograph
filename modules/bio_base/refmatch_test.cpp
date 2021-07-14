#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/refmatch.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/io/spiral_file_mem.h"

using namespace testing;
using namespace dna_testutil;

class refmatch_test : public Test {
 protected:
  void make_rmat() {
    spiral_file_create_mem c;

    refmatch_builder b(m_seqset, m_ref.get());
    b.build(c.create());

    spiral_file_open_mem o(c.close());
    m_rmat.emplace(m_seqset, m_ref.get(), o.open());
  }

  void use_reference(const std::vector<dna_sequence>& seqs) {
    m_ref = create_reference(seqs);
  }

  void use_reads(const std::vector<dna_sequence>& reads) {
    m_seqset_f = seqset_for_reads(reads);
    m_seqset = &m_seqset_f->get_seqset();
  }

  refmatch::entry get_rmat(const dna_sequence& seq) {
    seqset_range r = m_seqset->find(seq);
    CHECK(r.valid()) << seq;
    CHECK_EQ(r.begin() + 1, r.end()) << seq;
    uint64_t seqset_id = r.begin();
    CHECK_EQ(int(r.size()), int(m_seqset->entry_size(seqset_id))) << seq;

    return m_rmat->get(seqset_id);
  }

  boost::optional<refmatch> m_rmat;
  std::unique_ptr<reference> m_ref;
  std::unique_ptr<seqset_file> m_seqset_f;
  const seqset* m_seqset = nullptr;
};

void update_progress(const float& new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

TEST_F(refmatch_test, simple) {
  use_reference({tseq("abcdefghijklmno"),
                 tseq("ABCDEFGHIjklmnop") + tseq_rc("mnop"), tseq("01234056"),
                 tseq_rc("560789")});
  use_reads({tseq("abcde"), tseq_rc("bcde"), tseq("lmno"), tseq("mnop"),
             tseq("mnOP"), tseq("op") + tseq_rc("op"), tseq("0123"),
             tseq("3056"), tseq("60789")});
  make_rmat();

  auto entry = get_rmat(tseq("abcde"));
  EXPECT_TRUE(entry.has_fwd());
  EXPECT_FALSE(entry.has_rev());
  EXPECT_EQ(1, entry.matches());

  entry = get_rmat(tseq_rc("abcde"));
  EXPECT_FALSE(entry.has_fwd());
  EXPECT_TRUE(entry.has_rev());
  EXPECT_EQ(1, entry.matches());

  entry = get_rmat(tseq("lmno"));
  EXPECT_TRUE(entry.has_fwd());
  EXPECT_FALSE(entry.has_rev());
  EXPECT_EQ(2, entry.matches());

  entry = get_rmat(tseq("mnop"));
  EXPECT_TRUE(entry.has_fwd());
  EXPECT_TRUE(entry.has_rev());
  EXPECT_EQ(2, entry.matches());

  entry = get_rmat(tseq("op") + tseq_rc("op"));
  EXPECT_TRUE(entry.has_fwd());
  EXPECT_TRUE(entry.has_rev());
  EXPECT_EQ(2, entry.matches());

  entry = get_rmat(tseq("mnOP"));
  EXPECT_FALSE(entry.has_fwd());
  EXPECT_FALSE(entry.has_rev());
  EXPECT_EQ(0, entry.matches());

  entry = get_rmat(tseq_rc("mnOP"));
  EXPECT_FALSE(entry.has_fwd());
  EXPECT_FALSE(entry.has_rev());
  EXPECT_EQ(0, entry.matches());

  entry = get_rmat(tseq("0789"));
  EXPECT_FALSE(entry.has_fwd());
  EXPECT_TRUE(entry.has_rev());
  EXPECT_EQ(1, entry.matches());
}

TEST_F(refmatch_test, overflow) {
  std::vector<std::pair<dna_sequence, size_t /* num repeats */>> repeats = {
      {tseq("abcde"), refmatch::k_count_mask - 2},
      {tseq("fghij"), refmatch::k_count_mask - 1},
      {tseq("klmno"), refmatch::k_count_mask},
      {tseq("pqrst"), refmatch::k_count_mask + 1},
      {tseq("uvwxy"), refmatch::k_count_mask + 2}};

  std::vector<dna_sequence> ref_seqs;
  std::vector<dna_sequence> reads;
  for (const auto& r : repeats) {
    dna_sequence seq;
    reads.push_back(r.first);
    for (unsigned i = 0; i < r.second; ++i) {
      seq += r.first;
    }
    ref_seqs.push_back(seq);
  }

  use_reference(ref_seqs);
  use_reads(reads);
  make_rmat();

  for (const auto& r : repeats) {
    auto entry = get_rmat(r.first);
    EXPECT_TRUE(entry.has_fwd()) << r.first;
    EXPECT_FALSE(entry.has_rev()) << r.first;
    EXPECT_EQ(r.second, entry.matches()) << r.first;
  }
}

// Make sure reference locations are counted exactly once, even when
// split up into chunks.
class refmatch_chunk_test
    : public refmatch_test,
      public WithParamInterface<
          std::tuple<std::pair<int /* offset */, int /* stride */>,
                     int /* min chunk size */>> {
  void SetUp() override {
    m_offset = std::get<0>(GetParam()).first;
    m_stride = std::get<0>(GetParam()).second;
    m_chunk_size = std::get<1>(GetParam());

    m_orig_chunk_size = refmatch_builder::g_min_chunk_size;
    refmatch_builder::g_min_chunk_size = m_chunk_size;
  }

  void TearDown() override {
    refmatch_builder::g_min_chunk_size = m_orig_chunk_size;
  }

 protected:
  unsigned m_offset;
  unsigned m_stride;
  unsigned m_chunk_size;
  unsigned m_orig_chunk_size;
};

TEST_P(refmatch_chunk_test, chunking) {
  size_t counter = 0;
  std::vector<dna_sequence> ref_seqs;
  std::vector<dna_sequence> read_seqs;
  constexpr unsigned k_counter_length = 3;
  for (size_t chunk_length : std::initializer_list<size_t>{
           100, m_chunk_size - 1, m_chunk_size, m_chunk_size + 1,
           m_chunk_size * 2 - 1, m_chunk_size * 2, m_chunk_size * 2 + 1}) {
    dna_sequence seq;
    size_t counter_bases = k_counter_length * k_dna_test_sequence_length;
    while (seq.size() < chunk_length) {
      size_t this_counter = counter++;
      std::string cstr;
      cstr.push_back(' ');
      for (unsigned i = 1; i < k_counter_length; ++i) {
        cstr.push_back(' ' + 1 + (this_counter % 94));
        this_counter /= 94;
      }

      ASSERT_EQ(this_counter, 0) << "k_counter_length too small";

      auto to_add = tseq(cstr);
      ASSERT_EQ(to_add.size(), counter_bases);
      seq += to_add;
    }
    if (seq.size() > chunk_length) {
      seq = seq.subseq(0, chunk_length);
    }
    ref_seqs.push_back(seq);

    size_t read_length = counter_bases + 3;
    for (size_t i = m_offset; i + read_length <= chunk_length; i += m_stride) {
      read_seqs.push_back(seq.subseq(i, read_length));
    }
  }

  use_reference(ref_seqs);
  use_reads(read_seqs);
  make_rmat();

  for (const dna_sequence& orig_seq : read_seqs) {
    {
      dna_sequence seq = orig_seq;
      auto entry = get_rmat(seq);
      EXPECT_TRUE(entry.has_fwd()) << seq;
      EXPECT_FALSE(entry.has_rev()) << seq;
      EXPECT_EQ(1, entry.matches()) << seq;
    }
    {
      dna_sequence seq = orig_seq.rev_comp();
      auto entry = get_rmat(seq);
      EXPECT_FALSE(entry.has_fwd()) << seq;
      EXPECT_TRUE(entry.has_rev()) << seq;
      EXPECT_EQ(1, entry.matches()) << seq;
    }
  }
}

INSTANTIATE_TEST_CASE_P(
    refmatch_chunk_tests, refmatch_chunk_test,
    ::testing::Combine(::testing::Values(std::make_pair(0, 1),
                                         std::make_pair(0, 2),
                                         std::make_pair(1, 2)),
                       ::testing::Values(1, 10, 999, 25600)));

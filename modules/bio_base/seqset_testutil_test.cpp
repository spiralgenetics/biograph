#include "modules/bio_base/seqset_testutil.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/readmap.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace dna_testutil;

}  // namespace

// Tests a simple case of seqset_for_reads
TEST(seqset_testutil_test, simple_seqset) {
  dna_sequence seq("GCTACGC");
  // TAGC
  auto sq_file = seqset_for_reads({seq});
  EXPECT_TRUE(sq_file->get_seqset().find(seq).valid());
  EXPECT_TRUE(sq_file->get_seqset().find(seq.rev_comp()).valid());
  EXPECT_EQ(10, sq_file->get_seqset().size());
}

TEST(seqset_testutil_test, simple_seqset2) {
  dna_sequence seq(tseq("a"));
  auto sq_file = seqset_for_reads({seq});
  EXPECT_TRUE(sq_file->get_seqset().find(seq).valid());
  EXPECT_TRUE(sq_file->get_seqset().find(seq.rev_comp()).valid());
  EXPECT_EQ(18, sq_file->get_seqset().size());
}

// Tests that a seqset was successfully generated, and that all the
// included sequences are present and all the excluded sequences are
// not.

TEST(seqset_testutil_test, seqset_for_reads) {
  std::vector<dna_sequence> test_seqs;
  for (int i = -128; i < 128; i++) {
    test_seqs.push_back(dna_test_sequence(std::string({char(i)})));
  }
  std::vector<dna_sequence> include_seqs;
  std::vector<dna_sequence> exclude_seqs;

  for (unsigned i = 0; i < test_seqs.size(); i++) {
    dna_sequence seq = test_seqs[i];
    if (i & 2) {
      seq = seq.rev_comp();
    }
    if (i & 1) {
      include_seqs.push_back(seq);
    } else {
      exclude_seqs.push_back(seq);
    }
  }

  std::unique_ptr<seqset_file> sq_file = seqset_for_reads(include_seqs);
  const seqset& sq = sq_file->get_seqset();

  for (unsigned i = 0; i < test_seqs.size(); i++) {
    dna_sequence seq = test_seqs[i];
    dna_sequence rc_seq = seq.rev_comp();

    if (i & 1) {
      EXPECT_TRUE(sq.find(seq).valid());
      EXPECT_TRUE(sq.find(rc_seq).valid());
    } else {
      EXPECT_FALSE(sq.find(seq).valid());
      EXPECT_FALSE(sq.find(rc_seq).valid());
    }
  }
  EXPECT_EQ(893, sq.size());
}

TEST(seqset_testutil_test, seqset_for_flats) {
  std::vector<dna_sequence> include_seqs;
  for (int i = -128; i < 128; i++) {
    include_seqs.push_back(dna_test_sequence(std::string({char(i)})));
  }

  std::unique_ptr<seqset_file> sq_file = seqset_for_reads(include_seqs);
  std::unique_ptr<seqset_flat> sq_flat = seqset_flat_for_seqset(&sq_file->get_seqset());

  for (size_t i = 0; i < sq_file->get_seqset().size(); ++i) {
    EXPECT_EQ(sq_file->get_seqset().ctx_entry(i).sequence().as_string(),
              sq_flat->get(i).as_string());
  }
  EXPECT_EQ(sq_flat->size(), sq_file->get_seqset().size());
}

TEST(seqset_testutil_test, long_reads) {
  dna_sequence read_source = tseq("abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  ASSERT_GT(read_source.size(), 300);

  std::vector<dna_sequence> expected_seqs;
  for (size_t read_idx = 0; read_idx < read_source.size(); ++read_idx) {
    expected_seqs.push_back((read_source + read_source).subseq(read_idx, read_source.size()));
  }
  std::shared_ptr<seqset> ss = seqset_for_reads(expected_seqs);
  std::unique_ptr<readmap> rm = readmap_for_reads(ss, {}, expected_seqs);

  std::vector<dna_sequence> actual;
  for (uint64_t seqset_id = 0; seqset_id < ss->size(); ++seqset_id) {
    actual.push_back(ss->ctx_entry(seqset_id).sequence());
  }

  std::vector<dna_sequence> actual_reads;
  for (uint32_t read_id = 0; read_id < rm->size(); ++read_id) {
    actual_reads.push_back(
        ss->ctx_entry(rm->index_to_entry(read_id)).sequence(rm->get_readlength(read_id)));
  }

  for (size_t read_idx = 0; read_idx < read_source.size(); ++read_idx) {
    expected_seqs.push_back(
        (read_source + read_source).rev_comp().subseq(read_idx, read_source.size()));
  }
  std::sort(expected_seqs.begin(), expected_seqs.end());

  EXPECT_THAT(actual, testing::ElementsAreArray(expected_seqs));
  EXPECT_THAT(actual_reads, testing::ElementsAreArray(expected_seqs));
}

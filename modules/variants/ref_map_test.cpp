#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/variants/ref_map.h"

using namespace testing;
using namespace dna_testutil;
using namespace variants;

MATCHER_P2(AnchorIs, scaffold_id, offset, "") {
  if (!arg) {
    return false;
  }
  auto anchor = *arg;
  return anchor.rev_comp == false && anchor.pos.scaffold_id == scaffold_id &&
         anchor.pos.position == (unsigned long)offset;
}

MATCHER_P2(AnchorIsRc, scaffold_id, offset, "") {
  if (!arg) {
    return false;
  }
  auto anchor = *arg;
  return anchor.rev_comp == true && anchor.pos.scaffold_id == scaffold_id &&
         anchor.pos.position == (unsigned long)offset;
}

namespace variants {

class ref_map_test : public Test {
 protected:
  void make_rmap() {
    m_rmap.emplace(m_seqset, m_ref.get());
    m_rmap->build();
  }

  void use_reference(const std::vector<dna_sequence>& seqs) {
    m_ref = create_reference(seqs);
  }

  void use_reads(const std::vector<dna_sequence>& reads) {
    m_seqset_f = seqset_for_reads(reads);
    m_seqset = &m_seqset_f->get_seqset();
  }

  ref_map::entry get_rmap(const dna_sequence& seq) {
    seqset_range r = m_seqset->find(seq);
    CHECK(r.valid()) << seq;

    CHECK_EQ(r.begin() + 1, r.end()) << seq;

    m_anchor = m_rmap->get_unique_ref_anchor(r.begin());
    return m_rmap->get(r.begin());
  }

  boost::optional<ref_map> m_rmap;
  std::unique_ptr<reference> m_ref;
  std::unique_ptr<seqset_file> m_seqset_f;
  const seqset* m_seqset = nullptr;
  static constexpr size_t k_min_chunk_size = ref_map::k_min_chunk_size;
  boost::optional<ref_anchor> m_anchor;
};

constexpr size_t ref_map_test::k_min_chunk_size;

TEST_F(ref_map_test, simple) {
  use_reference({tseq("abcdefghijklmno"),
                 tseq("ABCDEFGHIjklmnop") + tseq_rc("mnop"), tseq("01234056"),
                 tseq_rc("560789")});
  use_reads({tseq("abcde"), tseq_rc("bcde"), tseq("lmno"), tseq("mnop"),
             tseq("mnOP"), tseq("op") + tseq_rc("op"), tseq("0123"),
             tseq("3056"), tseq("60789")});
  make_rmap();

  auto entry = get_rmap(tseq("abcde"));
  EXPECT_TRUE(entry.fwd_match());
  EXPECT_FALSE(entry.rev_match());
  EXPECT_EQ(1, entry.match_count());
  EXPECT_THAT(m_anchor, AnchorIs(0, 0));

  entry = get_rmap(tseq_rc("abcde"));
  EXPECT_FALSE(entry.fwd_match());
  EXPECT_TRUE(entry.rev_match());
  EXPECT_EQ(1, entry.match_count());
  EXPECT_THAT(m_anchor, AnchorIsRc(0, k_dna_test_sequence_length * 5));

  entry = get_rmap(tseq("bcde"));
  EXPECT_TRUE(entry.fwd_match());
  EXPECT_FALSE(entry.rev_match());
  EXPECT_EQ(1, entry.match_count());
  EXPECT_THAT(m_anchor, AnchorIs(0, k_dna_test_sequence_length));

  entry = get_rmap(tseq_rc("bcde"));
  EXPECT_FALSE(entry.fwd_match());
  EXPECT_TRUE(entry.rev_match());
  EXPECT_EQ(1, entry.match_count());
  EXPECT_THAT(m_anchor, AnchorIsRc(0, k_dna_test_sequence_length * 5));

  entry = get_rmap(tseq("lmno"));
  EXPECT_TRUE(entry.fwd_match());
  EXPECT_FALSE(entry.rev_match());
  EXPECT_EQ(2, entry.match_count());
  EXPECT_EQ(m_anchor, boost::none);

  entry = get_rmap(tseq("mnop"));
  EXPECT_TRUE(entry.fwd_match());
  EXPECT_TRUE(entry.rev_match());
  EXPECT_EQ(2, entry.match_count());
  EXPECT_EQ(m_anchor, boost::none);

  entry = get_rmap(tseq("op") + tseq_rc("op"));
  EXPECT_TRUE(entry.fwd_match());
  EXPECT_TRUE(entry.rev_match());
  EXPECT_EQ(2, entry.match_count());
  EXPECT_EQ(m_anchor, boost::none);

  entry = get_rmap(tseq("op"));
  EXPECT_TRUE(entry.fwd_match());
  EXPECT_TRUE(entry.rev_match());
  EXPECT_EQ(2, entry.match_count());
  EXPECT_EQ(m_anchor, boost::none);

  entry = get_rmap(tseq("mnOP"));
  EXPECT_FALSE(entry.fwd_match());
  EXPECT_FALSE(entry.rev_match());
  EXPECT_EQ(0, entry.match_count());
  EXPECT_EQ(m_anchor, boost::none);

  entry = get_rmap(tseq_rc("mnOP"));
  EXPECT_FALSE(entry.fwd_match());
  EXPECT_FALSE(entry.rev_match());
  EXPECT_EQ(0, entry.match_count());
  EXPECT_EQ(m_anchor, boost::none);

  entry = get_rmap(tseq("01"));
  EXPECT_TRUE(entry.fwd_match());
  EXPECT_FALSE(entry.rev_match());
  EXPECT_EQ(1, entry.match_count());
  EXPECT_THAT(m_anchor, AnchorIs(2, 0));

  entry = get_rmap(tseq("0789"));
  EXPECT_FALSE(entry.fwd_match());
  EXPECT_TRUE(entry.rev_match());
  EXPECT_EQ(1, entry.match_count());
  EXPECT_THAT(m_anchor, AnchorIsRc(3, 4 * k_dna_test_sequence_length));
}

TEST_F(ref_map_test, get_ref_slice) {
  use_reference({tseq("abcdefghijklmno"),
                 tseq("ABCDEFGHIjklmnop") + tseq_rc("mnop"), tseq("01234056"),
                 tseq_rc("560789")});

  use_reads({tseq("abcde")});
  make_rmap();
  ref_anchor a;
  a.rev_comp = false;
  a.pos.position = tseq("abcde").size();
  a.pos.scaffold_id = 0;
  EXPECT_EQ(m_rmap->get_ref_slice(a), tseq("fghijklmno"));
  a.pos.position = tseq("ABCDE").size();
  a.pos.scaffold_id = 1;
  a.rev_comp = true;
  EXPECT_EQ(m_rmap->get_ref_slice(a), tseq("ABCDE").rev_comp());
}

// Make sure reference locations are counted exactly once, even when
// split up into chunks.
class ref_map_chunk_test
    : public ref_map_test,
      public WithParamInterface<std::pair<int /* offset */, int /* stride */>> {
};

TEST_P(ref_map_chunk_test, chunking) {
  auto offset = std::get<0>(GetParam());
  auto stride = std::get<1>(GetParam());

  size_t counter = 0;
  constexpr int k_counter_length = 3;
  std::vector<dna_sequence> ref_seqs;
  std::vector<dna_sequence> read_seqs;
  for (size_t chunk_length : {100UL, k_min_chunk_size - 1, k_min_chunk_size,
                              k_min_chunk_size + 1, k_min_chunk_size * 2 - 1,
                              k_min_chunk_size * 2, k_min_chunk_size * 2 + 1}) {
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
    for (size_t i = offset; i + read_length <= chunk_length; i += stride) {
      read_seqs.push_back(seq.subseq(i, read_length));
    }
  }

  use_reference(ref_seqs);
  use_reads(read_seqs);
  make_rmap();

  for (const dna_sequence& orig_seq : read_seqs) {
    auto stride_counter = stride;

    dna_sequence seq = orig_seq;
    while (stride_counter) {
      auto entry = get_rmap(seq);
      EXPECT_TRUE(entry.fwd_match()) << seq;
      EXPECT_FALSE(entry.rev_match()) << seq;
      EXPECT_EQ(1, entry.match_count()) << seq;

      --stride_counter;
      seq = seq.subseq(0, seq.size() - 1);
    }

    seq = orig_seq.rev_comp();
    while (stride_counter) {
      auto entry = get_rmap(seq.rev_comp());
      EXPECT_FALSE(entry.fwd_match()) << seq;
      EXPECT_TRUE(entry.rev_match()) << seq;
      EXPECT_EQ(1, entry.match_count()) << seq;

      --stride_counter;
      seq = seq.subseq(0, seq.size() - 1);
    }
  }
}

INSTANTIATE_TEST_CASE_P(ref_map_chunk_tests, ref_map_chunk_test,
                        ::testing::Values(std::make_pair(0, 1),
                                          std::make_pair(0, 2),
                                          std::make_pair(1, 2)));

}  // namespace variants

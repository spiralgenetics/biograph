#include "modules/variants/scaffold.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

// Copy extents into an array of pairs for easier matching
std::vector<std::pair<aoffset_t, dna_sequence>> get_extents(const scaffold& s) {
  std::vector<std::pair<aoffset_t, dna_sequence>> out;
  for (const auto& e : s.extents()) {
    out.push_back(std::make_pair(e.offset, dna_sequence(e.sequence.begin(), e.sequence.end())));
  }
  return out;
}

TEST(scaffold_test, simple) {
  scaffold s;
  EXPECT_FALSE(s.is_simple());
  EXPECT_TRUE(s.empty());
  EXPECT_EQ(0, s.end_pos());
  s.add(0, tseq("abcde"));
  EXPECT_TRUE(s.is_simple());
  EXPECT_EQ(s.end_pos(), tseq("abcde").size());
  EXPECT_EQ(s.get_simple(), tseq("abcde"));
  EXPECT_FALSE(s.empty());
}

TEST(scaffold_test, two_regions) {
  scaffold s;
  s.add(5, tseq("abcde"));
  s.add(100, tseq("fghi"));

  EXPECT_EQ(100 + tseq("fghi").size(), s.end_pos());
  EXPECT_FALSE(s.is_simple());
  EXPECT_FALSE(s.empty());
}

TEST(scaffold_test, empty_subscaffold) {
  scaffold s;
  s.add(5, tseq("abcde"));
  s.add(100, tseq("fghi"));

  scaffold sub = s.subscaffold(80, 1);
  EXPECT_TRUE(sub.empty());
  EXPECT_FALSE(sub.is_simple());
}

TEST(scaffold_test, offset_subscaffold) {
  scaffold s;
  s.add(5, tseq("abcde"));
  s.add(100, tseq("fghi"));

  scaffold sub = s.subscaffold(1, 100 + tseq("fghi").size() - 1);
  EXPECT_EQ(sub.end_pos(), 100 + tseq("fghi").size() - 1);
  EXPECT_THAT(get_extents(sub), ElementsAre(Pair(4, tseq("abcde")), Pair(99, tseq("fghi"))));
  EXPECT_FALSE(sub.empty());
}

TEST(scaffold_test, partial_subscaffold) {
  scaffold s;
  s.add(5, tseq("abcde"));
  s.add(100, tseq("fghi"));

  scaffold sub = s.subscaffold(10, 98);
  EXPECT_EQ(sub.end_pos(), 98);
  EXPECT_THAT(get_extents(sub),
              ElementsAre(Pair(0, tseq("abcde").subseq(10 - 5, tseq("abcde").size() - (10 - 5))),
                          Pair(100 - 10, tseq("fghi").subseq(0, 8))));
}

TEST(scaffold_test, subscaffold_inside) {
  scaffold s;
  s.add(5, tseq("abcde"));
  s.add(100, tseq("fghi"));

  scaffold sub = s.subscaffold(5 + tseq("abcde").size(), 100 - (5 + tseq("abcde").size()));
  EXPECT_TRUE(sub.empty());
  EXPECT_EQ(100 - (5 + tseq("abcde").size()), sub.end_pos());
  EXPECT_FALSE(sub.is_simple());
}

TEST(scaffold_test, subscaffold_outside) {
  scaffold s;
  s.add(5, tseq("abcde"));
  s.add(100, tseq("fghi"));

  scaffold sub = s.subscaffold(5, 100 + tseq("fghi").size() - 5);
  EXPECT_EQ(100 + tseq("fghi").size() - 5, sub.end_pos());
  EXPECT_FALSE(sub.empty());
  EXPECT_THAT(get_extents(sub), ElementsAre(Pair(0, tseq("abcde")), Pair(100 - 5, tseq("fghi"))));
  EXPECT_FALSE(sub.is_simple());
}

TEST(scaffold_test, as_string) {
  scaffold s;
  s.add(0, tseq("a"));
  s.add(tseq("a").size() + 5, tseq("b"));
  s.set_end_pos(tseq("a").size() + 5 + tseq("b").size() + 3);

  EXPECT_EQ(s.as_string(), tseq("a").as_string() + "NNNNN" + tseq("b").as_string() + "NNN");
}

TEST(scaffold_test, iterator) {
  scaffold s;
  s.add(5, dna_sequence("GA"));
  s.add(10, dna_sequence("T"));

  auto cur = s.begin();
  auto end = s.end();
  ASSERT_FALSE(cur == end);

  EXPECT_EQ(5, cur.offset());
  EXPECT_EQ(dna_base('G'), *cur);
  EXPECT_EQ(dna_base('C'), cur->complement());

  ++cur;
  ASSERT_FALSE(cur == end);

  EXPECT_EQ(6, cur.offset());
  EXPECT_EQ(dna_base('A'), *cur);
  EXPECT_EQ(dna_base('T'), cur->complement());

  ++cur;
  ASSERT_FALSE(cur == end);

  EXPECT_EQ(10, cur.offset());
  EXPECT_EQ(dna_base('T'), *cur);
  EXPECT_EQ(dna_base('A'), cur->complement());

  ++cur;
  EXPECT_TRUE(cur == end);
  EXPECT_EQ(11, cur.offset());
}

TEST(scaffold_test, skip) {
  scaffold s;
  s.add(0, dna_sequence("C"));
  s.add(5, dna_sequence("GA"));
  s.add(10, dna_sequence("T"));
  s.set_end_pos(15);

  auto it = s.begin();
  EXPECT_EQ(it.offset(), 0);
  it.skip_to(4, "scaffold_test");
  EXPECT_EQ(it.offset(), 5);
  EXPECT_EQ(dna_base('G'), *it);
  EXPECT_TRUE(it.first_in_extent());

  it = s.begin();
  it.skip_to(5, "scaffold_test");
  EXPECT_EQ(it.offset(), 5);
  EXPECT_EQ(dna_base('G'), *it);
  EXPECT_TRUE(it.first_in_extent());

  it = s.begin();
  it.skip_to(6, "scaffold_test");
  EXPECT_EQ(it.offset(), 6);
  EXPECT_EQ(dna_base('A'), *it);
  EXPECT_FALSE(it.first_in_extent());

  it = s.begin();
  it.skip_to(7, "scaffold_test");
  EXPECT_EQ(it.offset(), 10);
  EXPECT_EQ(dna_base('T'), *it);
  EXPECT_TRUE(it.first_in_extent());

  it = s.begin();
  it.skip_to(10, "scaffold_test");
  EXPECT_EQ(it.offset(), 10);
  EXPECT_EQ(dna_base('T'), *it);
  EXPECT_TRUE(it.first_in_extent());

  it = s.begin();
  it.skip_to(11, "scaffold_test");
  EXPECT_EQ(it.offset(), 15);
  EXPECT_TRUE(it == s.end());
}

}  // namespace variants

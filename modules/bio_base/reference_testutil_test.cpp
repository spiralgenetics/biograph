#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/dna_testutil.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;

MATCHER_P2(ScaffoldIs, name, size, "") {
  return arg.name == std::to_string(name) && arg.size == size;
}

MATCHER_P4(ExtentIs, scaffold_name, size, offset, flat, "") {
  return arg.scaffold_name == size_t(scaffold_name) && arg.size == size &&
         arg.offset == size_t(offset) && arg.flat == flat;
}

TEST(reference_testutil_test, create_flat_ref) {
  std::vector<dna_sequence> seqs = {
      dna_test_sequence("x") + dna_test_sequence("y"),
      dna_test_sequence("x") + dna_sequence("A") + dna_test_sequence("y"),
      dna_test_sequence("x") + dna_sequence("AA") + dna_test_sequence("y")};

  const size_t size1 = seqs[0].size();
  const size_t size2 = seqs[1].size();
  const size_t size3 = seqs[2].size();

  // Apparently flat offsets start at 1? TODO(nils): Investigate
  // further and document more clearly.
  const size_t offset1 = 1;
  const size_t offset2 = offset1 + size1;
  const size_t offset3 = offset2 + size2;

  std::unique_ptr<flat_ref> ref = create_flat_ref(seqs);

  const flat_ref::index_t& index = ref->get_index();

  EXPECT_THAT(index.scaffolds, ElementsAre(ScaffoldIs(0, size1),  //
                                           ScaffoldIs(1, size2),  //
                                           ScaffoldIs(2, size3)));

  EXPECT_THAT(index.extents, ElementsAre(ExtentIs(0, size1, 0, offset1),
                                         ExtentIs(1, size2, 0, offset2),
                                         ExtentIs(2, size3, 0, offset3)));

  auto it1 = ref->get_dna(offset1);
  auto it2 = ref->get_dna(offset2);
  auto it3 = ref->get_dna(offset3);

  dna_sequence result1(it1, it1 + size1);
  dna_sequence result2(it2, it2 + size2);
  dna_sequence result3(it3, it3 + size3);

  EXPECT_EQ(seqs[0], result1);
  EXPECT_EQ(seqs[1], result2);
  EXPECT_EQ(seqs[2], result3);
}

TEST(reference_testutil_test, create_reference) {
  std::vector<dna_sequence> seqs = {
      dna_test_sequence("x") + dna_test_sequence("y"),
      dna_test_sequence("x") + dna_sequence("A") + dna_test_sequence("y"),
      dna_test_sequence("x") + dna_sequence("AA") + dna_test_sequence("y")};

  std::unique_ptr<reference> full_ref = create_reference(seqs);
  const flat_ref& ref = full_ref->get_flat_ref();
  auto bwt = full_ref->get_bwt();
  const flat_ref::index_t& index = ref.get_index();

  std::vector<dna_slice> ref_seqs;
  size_t offset = 0;
  std::map<dna_sequence, size_t> expected_at;
  std::set<size_t> all_expected_matches;
  for (const flat_ref::extent_t& extent : index.extents) {
    auto data_start = ref.get_dna(extent.flat);
    dna_slice slice = dna_slice(data_start, data_start + extent.size);
    ref_seqs.push_back(slice);
    if (slice.size() > 1) {
      expected_at[dna_sequence(slice.begin(), slice.end())] = offset;
      all_expected_matches.insert(offset);
    }
    offset += slice.size();
  }

  EXPECT_THAT(ref_seqs,
              UnorderedElementsAre(seqs[0], seqs[1], seqs[2],
                                   // The testutil also includes one
                                   // of each base to avoid bwt
                                   // building errors due to lack of
                                   // at least 1 of each base.
                                   dna_sequence("A"), dna_sequence("C"),
                                   dna_sequence("G"), dna_sequence("T")));

  for (const auto& e_at : expected_at) {
    const dna_sequence& seq = e_at.first;
    size_t offset = e_at.second;

    if (seq.size() == 1) {
      continue;
    }

    auto found = bwt.find(seq);
    EXPECT_TRUE(found.valid());
    if (!found.valid()) {
      continue;
    }

    EXPECT_EQ(1, found.matches()) << seq;
    EXPECT_EQ(offset, found.get_match(0)) << seq;
  }

  // All the reference sequences start with "x".
  auto found = bwt.find(dna_test_sequence("x"));
  ASSERT_TRUE(found.valid());
  EXPECT_EQ(all_expected_matches.size(), found.matches());
  std::set<size_t> all_actual_matches;
  for (size_t i = 0; i != found.matches(); ++i) {
    all_actual_matches.insert(found.get_match(i));
  }
  EXPECT_THAT(all_actual_matches, ContainerEq(all_expected_matches));

  auto not_found = bwt.find(dna_test_sequence("y") + dna_test_sequence("x"));
  EXPECT_FALSE(not_found.valid());
}

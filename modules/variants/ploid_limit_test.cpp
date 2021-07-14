#include "modules/variants/ploid_limit.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class ploid_test : public assemble_test {
 public:
  void add_ref(int id, aoffset_t left, aoffset_t right, acost_t score) {
    add(id, left, right, score, true /* ref matches */);
  }

  void add(int id, aoffset_t left, aoffset_t right, acost_t score, bool ref_matches = false,
           aoffset_t left_anchor_len = 0, aoffset_t right_anchor_len = 0) {
    assembly a;
    a.assembly_id = id;
    a.left_offset = left;
    a.right_offset = right;
    CHECK_LE(left + left_anchor_len + right_anchor_len, right);
    a.score = score;
    a.matches_reference = ref_matches;
    if (ref_matches) {
      dna_slice seq = m_scaffold.subscaffold(left, right - left).get_simple();
      a.seq = dna_sequence(seq.begin(), seq.end());
    } else {
      std::uniform_int_distribution<aoffset_t> var_len(
          0, (right - left - right_anchor_len - left_anchor_len) * 2 / 10);
      std::uniform_int_distribution<int> is_zero(0, 1);
      aoffset_t seq_size = 0;
      if (!is_zero(*m_rand_source)) {
        seq_size = var_len(*m_rand_source) * 10;
      }
      aligned_var v;
      v.left_offset = a.left_offset + left_anchor_len;
      v.right_offset = a.right_offset - right_anchor_len;
      while (int(v.seq.size()) < seq_size) {
        v.seq += tseq("abcdefghijklmnopqrstuvwxyz");
      }
      v.seq = v.seq.subseq(0, seq_size);
      a.aligned_variants.push_back(v);

      if (left_anchor_len) {
        dna_slice left_anchor = m_scaffold.subscaffold(left, left_anchor_len).get_simple();
        a.seq += dna_sequence(left_anchor.begin(), left_anchor.end());
      }
      a.seq += v.seq;
      if (right_anchor_len) {
        dna_slice right_anchor =
            m_scaffold.subscaffold(right - right_anchor_len, right_anchor_len).get_simple();
        a.seq += dna_sequence(right_anchor.begin(), right_anchor.end());
      }
    }
    m_limiter->add(make_unique<assembly>(a));
  }

  void flush() {
    m_limiter.reset();
    expect_sorted(assembly::left_offset_less_than);
  }

  void SetUp() {
    assemble_options m_options;
    m_limiter = make_unique<ploid_limiter>(m_options, test_output());
    m_limiter->set_max_ploids(2);

    dna_sequence test_seq;
    while (test_seq.size() < 10000) {
      test_seq += tseq("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    }
    m_scaffold.add(0, test_seq);
    m_options.scaffold = &m_scaffold;

    std::random_device rand_dev;
    auto seed = rand_dev();
    std::cout << "Random seed is " << seed << "\n";
    m_rand_source.emplace(seed);
  }

 protected:
  std::unique_ptr<ploid_limiter> m_limiter;
  boost::optional<std::mt19937> m_rand_source;
};

MATCHER_P(AssemblyIdsAreArray, expected_ids,
          std::string(negation ? "doesn't have assembly ids " : "has assembly ids ") +
              PrintToString(expected_ids)) {
  std::vector<size_t> ids;
  ids.push_back(arg.assembly_id);
  for (const auto& id : arg.merged_assembly_ids) {
    ids.push_back(id);
  }

  return ExplainMatchResult(UnorderedElementsAreArray(expected_ids), ids, result_listener);
}

template <typename... Arg>
Matcher<const assembly&> AssemblyIdsAre(Arg... args) {
  return AssemblyIdsAreArray(std::vector<size_t>{size_t(args)...});
}

TEST_F(ploid_test, single_at_once) {
  add(1, 10, 20, 0);
  add(2, 30, 40, 0);

  // First one should come out once we get the second one, since we don't have to test for merging.
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(1)));

  flush();
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(1), AssemblyIdIs(2)));
}

TEST_F(ploid_test, single_ref_at_once) {
  add_ref(1, 10, 20, 0);
  add_ref(2, 30, 40, 0);

  // First ref assembly should come out without waiting for second
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(1)));

  flush();
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(1), AssemblyIdIs(2)));
}

TEST_F(ploid_test, two_at_time) {
  add(1, 10, 30, 0);
  add(2, 20, 40, 0);

  // Nothing should come out yet, since we don't know if anything will
  // have a better score.
  EXPECT_THAT(m_assemblies, IsEmpty());

  // Advancing past the end should output the assemblies.
  add_ref(3, 50, 60, 0);
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(1), AssemblyIdIs(2)));

  flush();
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(1), AssemblyIdIs(2), AssemblyIdIs(3)));
}

TEST_F(ploid_test, two_ref_at_time) {
  add_ref(1, 10, 30, 0);
  add_ref(2, 20, 40, 0);

  // Nothing should come out yet, since we don't know if anything will
  // have a better score.
  EXPECT_THAT(m_assemblies, IsEmpty());

  flush();
}

TEST_F(ploid_test, three_at_a_time) {
  add(1, 10, 40, 0);
  add(2, 20, 50, 100);
  add(3, 30, 60, 4);

  // Nothing should come out yet, since we don't know if anything will
  // have a better score.
  EXPECT_THAT(m_assemblies, IsEmpty());

  flush();
  // Only the two higher scoring ones should come out
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(2), AssemblyIdIs(3)));
}

TEST_F(ploid_test, three_ref_at_a_time) {
  add_ref(1, 10, 40, 0);
  add_ref(2, 20, 50, 100);
  add_ref(3, 30, 60, 4);

  // Nothing should come out yet, since we don't know if anything will
  // have a better score.
  EXPECT_THAT(m_assemblies, IsEmpty());

  add_ref(4, 50, 60, 0);
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(1), AssemblyIdIs(2)));

  flush();

  // Ref matches don't conflict with each other.
  EXPECT_THAT(m_assemblies,
              ElementsAre(AssemblyIdIs(1), AssemblyIdIs(2), AssemblyIdIs(3), AssemblyIdIs(4)));
}

TEST_F(ploid_test, three_but_not_at_a_time) {
  // Sacrifical entry to make sure the whole thing goes through the deploding
  // process
  add_ref(1, 0, 1000, 0);
  // Now, three at a time, but not all at once.
  add(2, 50, 500, 100);
  add(3, 100, 900, 100);
  add(4, 510, 950, 100);

  // Nothing should come out yet, since we don't know if anything will
  // have a better score.
  EXPECT_THAT(m_assemblies, IsEmpty());

  // This ref still overlaps, nothing should get output until we're done with variants that might
  // overlap anything.
  add_ref(5, 940, 1100, 0);
  EXPECT_THAT(m_assemblies, IsEmpty());

  add_ref(6, 1000, 1200, 0);
  // Now we're past all the variants, so they should get output.
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(2), AssemblyIdIs(3), AssemblyIdIs(4)));

  flush();
  // Ref matches don't conflict with each other.
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(2), AssemblyIdIs(3), AssemblyIdIs(4),
                                        AssemblyIdIs(5), AssemblyIdIs(6)));
}

TEST_F(ploid_test, three_with_one_better) {
  // Sacrifical entry to make sure the whole thing goes through the deploding
  // process
  add_ref(1, 0, 1000, 0);
  // Now, three at a time, but not all at once.
  add(2, 50, 500, 100);
  add(3, 100, 900, 200);
  add(4, 510, 950, 100);

  // This should make 4 go away, but no others.
  add(5, 800, 1000, 500);

  // Nothing should come out yet, since we don't know if anything will
  // have a better score.
  EXPECT_THAT(m_assemblies, IsEmpty());

  flush();

  // Ref matches don't conflict with each other.
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdIs(2), AssemblyIdIs(3), AssemblyIdIs(5)));
}

TEST_F(ploid_test, merge_with_ref) {
  add(1, 0, 50, 0);
  add_ref(2, 0, 20, 1);
  add(3, 10, 50, 200, false /* not ref */, 10, 10);
  add_ref(4, 40, 60, 2);

  EXPECT_THAT(m_assemblies, IsEmpty());
  add_ref(5, 50, 70, 5);

  ASSERT_THAT(m_assemblies, ElementsAre(AssemblyIdsAre(2, 3, 4), AssemblyIdIs(1)));
  const auto& a = m_assemblies.front();
  EXPECT_EQ(a.left_offset, 0);
  EXPECT_EQ(a.right_offset, 60);
  ASSERT_EQ(a.aligned_variants.size(), 1);
  const auto& v = a.aligned_variants.front();
  EXPECT_EQ(v.left_offset, 20);
  EXPECT_EQ(v.right_offset, 40);

  flush();
  EXPECT_THAT(m_assemblies, ElementsAre(AssemblyIdsAre(2, 3, 4), AssemblyIdIs(1), AssemblyIdIs(5)));
}

}  // namespace variants

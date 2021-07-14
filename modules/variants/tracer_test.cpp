#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/ref_map.h"
#include "modules/variants/reversable_tracer.h"

#include <fstream>

using namespace testing;
using namespace dna_testutil;

namespace variants {

MATCHER_P3(TrAssemblyIs, left_anchor, expected_seq, right_anchor, "") {
  return aoffset_t(arg.seq.size()) >
             (arg.left_anchor_len + arg.right_anchor_len) &&
         arg.left_anchor_len == int(left_anchor.size()) &&
         arg.right_anchor_len == int(right_anchor.size()) &&
         arg.seq == (left_anchor + expected_seq + right_anchor);
}

MATCHER_P(RefTrAssemblyIs, ref_seq, "") {
  return (arg.right_offset - arg.left_offset) == int(ref_seq.size()) &&
         arg.left_anchor_len == 0 && arg.right_anchor_len == 0 &&
         arg.matches_reference && arg.seq == ref_seq;
}

class tracer_test : public TestWithParam<std::pair<bool /* rev_comp */, bool /* use pop tracer */>> {
 protected:
  void use_ref_parts(std::vector<std::pair<aoffset_t, dna_sequence>> parts) {
    CHECK(m_scaffold.empty());
    aoffset_t last_pos = 0;
    for (const auto& i : parts) {
      EXPECT_GE(i.first, last_pos);
      m_scaffold.add(i.first, i.second);
      last_pos = i.first + i.second.size();
    }

    for (const auto& i : parts) {
      CHECK_EQ(i.second, get_ref_part_seq(i.first, i.second.size())) << " at "
                                                                     << i.first;
    }
  }

  dna_sequence get_ref_part_seq(aoffset_t offset, int len) {
    std::cout << "Getting part seq for " << offset << ", len " << len << "\n";
    scaffold sub = m_scaffold.subscaffold(offset, len);
    dna_slice slice;
    if (sub.is_simple()) {
      slice = sub.get_simple();
    } else {
      EXPECT_TRUE(sub.is_simple()) << sub << " from " << offset << " to "
                                   << len;
    }
    return dna_sequence(slice.begin(), slice.end());
  }

  void use_reads(std::vector<std::pair<dna_sequence, dna_sequence>> pairs,
                 std::vector<dna_sequence> reads = {}) {
    std::vector<dna_sequence> all_reads;
    for (const auto& pair : pairs) {
      all_reads.push_back(pair.first);
      all_reads.push_back(pair.second);
    }
    for (const auto& read : reads) {
      all_reads.push_back(read);
    }
    m_seqset = seqset_for_reads(all_reads);
    m_readmap = readmap_for_reads(m_seqset, pairs, reads);

    m_rmap.emplace(m_seqset.get(), m_ref.get());
    m_rmap->build();

    m_assemblies.clear();
    m_opts.seqset = m_seqset.get();
    m_opts.readmap = m_readmap.get();
    m_opts.ref = m_ref.get();
    m_opts.rmap = &m_rmap.get();
    m_opts.use_bidir_tracer = m_use_bidir_tracer;
    m_opts.pop_trace_anchor_drop = false;
    g_trace_all_assemblies = true;
    // TODO(nils): Remove this option and rework the test not to need it.
    m_opts.trace_reference_assemblies = true;
    m_output_f = make_unique<assemble_lambda_output>(
        [&](assembly_ptr a) {
          static std::mutex mu;
          std::lock_guard<std::mutex> l(mu);
          std::cout << "Got assembly: " << *a << "\n";
          m_assemblies.push_back(*a);
          if (m_rev_comp) {
            m_left_pair_matches.push_back(a->right_pair_matches.size());
            m_right_pair_matches.push_back(a->left_pair_matches.size());
          } else {
            m_right_pair_matches.push_back(a->right_pair_matches.size());
            m_left_pair_matches.push_back(a->left_pair_matches.size());
          }

          ASSERT_GT(a->seq.size(), a->left_anchor_len);
          ASSERT_GT(a->seq.size(), a->right_anchor_len);
          if (a->left_anchor_len) {
            EXPECT_EQ(a->seq.subseq(0, a->left_anchor_len),
                      get_ref_part_seq(a->left_offset, a->left_anchor_len))
                << *a;
          }
          if (a->right_anchor_len) {
            EXPECT_EQ(
                a->seq.subseq(aoffset_t(a->seq.size()) - a->right_anchor_len,
                              a->right_anchor_len),
                get_ref_part_seq(a->right_offset - a->right_anchor_len,
                                 a->right_anchor_len))
                << *a;
          }
        },
        "raw_assemblies");
    m_opts.scaffold = &m_scaffold;

    m_trace.emplace(m_rev_comp, m_opts);
    assemble_stats st = m_trace->assemble(m_output_f.get());
    std::cout << "Stats: " << st << "\n";

    if (m_opts.debug_paths) {
      std::ofstream of("/tmp/path-debug.dot");
      m_trace->output_path_debug_dot(of);
      std::cout << "Outputted path debug to /tmp/path-debug.dot\n";
    }
    m_trace.reset();

    if (m_rev_comp) {
      std::reverse(m_assemblies.begin(), m_assemblies.end());
      std::reverse(m_right_pair_matches.begin(), m_right_pair_matches.end());
      std::reverse(m_left_pair_matches.begin(), m_left_pair_matches.end());
      std::swap(m_left_pair_matches, m_right_pair_matches);
    }
  }

  void SetUp() override {
    if (!m_ref) {
      m_ref = create_reference({tseq("abcdefghijklmnopqrstuvwxyz"),
                                tseq("ABCDEFGHIJKLM_12345_NOPQRSTUVWXYZ"),
                                tseq("0123456789")});
    }

    m_opts.min_overlap = k_dna_test_sequence_length * 2;
    std::tie(m_rev_comp, m_use_bidir_tracer) = GetParam();
  }

  assemble_options m_opts;
  pipeline_step_t m_output_f;
  bool m_rev_comp;
  bool m_use_bidir_tracer;
  scaffold m_scaffold;
  boost::optional<ref_map> m_rmap;
  static std::unique_ptr<reference> m_ref;
  std::shared_ptr<seqset> m_seqset;
  std::unique_ptr<readmap> m_readmap;
  std::vector<assembly> m_assemblies;
  std::vector<int> m_left_pair_matches, m_right_pair_matches;
  boost::optional<reversable_tracer> m_trace;
};

std::unique_ptr<reference> tracer_test::m_ref;

TEST_P(tracer_test, only_reference) {
  use_ref_parts({{5, tseq("abcdefghijklmnopqrst")}});
  use_reads({{tseq("abcde"), tseq_rc("ijklm")},
             {tseq("cdefg"), tseq_rc("klmno")},
             {tseq("efghi"), tseq_rc("mnopq")},
             {tseq("ghijk"), tseq_rc("opqrs")}});
  EXPECT_THAT(m_assemblies, ElementsAre(RefTrAssemblyIs(tseq("abcdefg")),
                                        RefTrAssemblyIs(tseq("cdefghi")),
                                        RefTrAssemblyIs(tseq("efghijk")),
                                        RefTrAssemblyIs(tseq("ghijklm")),
                                        RefTrAssemblyIs(tseq("ijklmno")),
                                        RefTrAssemblyIs(tseq("klmnopq")),
                                        RefTrAssemblyIs(tseq("mnopqrs"))));

  EXPECT_THAT(m_right_pair_matches, ElementsAre(1, 1, 1, 0, 0, 0, 0));
  EXPECT_THAT(m_left_pair_matches, ElementsAre(0, 0, 0, 0, 1, 1, 1));
}

TEST_P(tracer_test, simple_variant) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvw")}});
  use_reads({{tseq("bcdef"), tseq_rc("klmno")},
             {tseq("efgHi"), tseq_rc("nopqr")},
             {tseq("gHijk"), tseq_rc("qrstu")},
             {tseq("Hijkl"), tseq_rc("stuvw")}});
  EXPECT_THAT(
      m_assemblies,
      ElementsAre(TrAssemblyIs(tseq("bcdef"), tseq("gHij"), tseq("klmno")),
                  RefTrAssemblyIs(tseq("klmnopqr")),
                  RefTrAssemblyIs(tseq("nopqrstu")),
                  RefTrAssemblyIs(tseq("qrstuvw"))));
  EXPECT_THAT(m_right_pair_matches, ElementsAre(1, 0, 0, 0));
  EXPECT_THAT(m_left_pair_matches, ElementsAre(0, 0, 0, 1));
}

TEST_P(tracer_test, compound_variant) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvw")}});
  use_reads({// Ref case:
             {tseq("bcdef"), tseq_rc("ghijk")},
             {tseq("defgh"), tseq_rc("ijklm")},
             {tseq("cdefg"), tseq_rc("lmnop")},

             // Non-ref case:
             {tseq("efgHi"), tseq_rc("nopqr")},
             {tseq("gHijk"), tseq_rc("qrstu")},
             {tseq("Hijkl"), tseq_rc("stuvw")}});
  EXPECT_THAT(
      m_assemblies,
      UnorderedElementsAre(
          // Ref case:
          RefTrAssemblyIs(tseq("bcdefg")), RefTrAssemblyIs(tseq("cdefgh")),
          RefTrAssemblyIs(tseq("defghijk")), RefTrAssemblyIs(tseq("ghijklm")),
          RefTrAssemblyIs(tseq("ijklmnop")), RefTrAssemblyIs(tseq("lmnopqr")),
          RefTrAssemblyIs(tseq("nopqrstu")), RefTrAssemblyIs(tseq("qrstuvw")),

          // Non-ref case:
          TrAssemblyIs(tseq("cdefg"), tseq("H"), tseq("ijklm"))));
}

TEST_P(tracer_test, cross_ref_sections) {
  use_ref_parts(
      {{0, tseq("abcdef")},
       {tseq("abcdef").size() + tseq("gh").size(), tseq("ijklmnopqr")}});
  use_reads({{tseq("abcde"), tseq_rc("ghijk")},
             {tseq("cdefg"), tseq_rc("ijklm")},
             {tseq("efghi"), tseq_rc("klmno")}});
  EXPECT_THAT(m_assemblies, ElementsAre(TrAssemblyIs(tseq("abcde"), tseq("fgh"),
                                                     tseq("ijklm")),
                                        RefTrAssemblyIs(tseq("ijklmno"))));

  // These pairs aren't far enough apart to get counted
  EXPECT_THAT(m_right_pair_matches, ElementsAre(0, 0));
  EXPECT_THAT(m_left_pair_matches, ElementsAre(0, 0));
}

TEST_P(tracer_test, cross_ref_delete) {
  use_ref_parts({{0, tseq("abcdefg")}, {100, tseq("hijklmno")}});
  use_reads({{tseq("abcde"), tseq_rc("cdefg")},
             {tseq("efghi"), tseq_rc("ghijk")},
             {tseq("ijklm"), tseq_rc("klmno")}});
  EXPECT_THAT(m_assemblies,
              ElementsAre(RefTrAssemblyIs(tseq("abcdefg")),
                          TrAssemblyIs(tseq("cdefg"), tseq("h"), tseq("ijklm")),
                          RefTrAssemblyIs(tseq("ijklmno"))));

  // These pairs aren't far enough apart to get counted
  EXPECT_THAT(m_right_pair_matches, ElementsAre(0, 0, 0));
  EXPECT_THAT(m_left_pair_matches, ElementsAre(0, 0, 0));
}

TEST_P(tracer_test, dead_end) {
  m_opts.min_anchor_drop_overlap = tseq("j").size();

  m_opts.min_overlap = tseq("abcd").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({}, reads_for_seq(tseq("abcde") + dna_G + tseq("f") + dna_G +
                                  tseq("ghi") + dna_G + tseq("j"),
                              tseq("abcde").size(),
                              tseq("a").size()  // Generate a 5-tseq read every
                                                // 1-tseq.
                              ));

  if (m_rev_comp) {
    // This dead end is only reachable going forward.
    EXPECT_THAT(m_assemblies, IsEmpty());
  } else {
    EXPECT_THAT(m_assemblies, ElementsAre(TrAssemblyIs(
                                  tseq("abcde"), dna_G + tseq("f") + dna_G +
                                                     tseq("ghi") + dna_G,
                                  tseq("j"))));
  }
}

TEST_P(tracer_test, dead_end2) {
  m_opts.min_anchor_drop_overlap = tseq("j").size() + 1;

  m_opts.min_overlap = tseq("abcd").size();
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({}, reads_for_seq(tseq("abcde") + dna_G + tseq("f") + dna_G +
                                  tseq("ghi") + dna_G + tseq("j"),
                              tseq("abcde").size(),
                              tseq("a").size()  // Generate a 5-tseq read every
                                                // 1-tseq.
                              ));

  if (m_rev_comp) {
    // This dead end is only reachable going forward.
    EXPECT_THAT(m_assemblies, IsEmpty());
  } else {
    // min_anchor_drop_overlap is too long to match the "j".
    EXPECT_THAT(m_assemblies,
                ElementsAre(TrAssemblyIs(
                    tseq("abcde"), dna_G + tseq("f") + dna_G, tseq("ghi"))));
  }
}

// This tests a dead end that's a larger drop than would be found with
// min_anchor_drop_overlap; the anchor is of length min_overlap.
TEST_P(tracer_test, dead_end_big_drop) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({{tseq("bcdef"), tseq("efgH") + dna_G + tseq("y")},
             {tseq("defgH") + dna_G, tseq("fgH") + dna_G + tseq("yz")}});
  if (m_rev_comp) {
    // This dead end is only reachable going forward.
    EXPECT_THAT(m_assemblies, IsEmpty());
  } else {
    EXPECT_THAT(m_assemblies,
                ElementsAre(TrAssemblyIs(tseq("bcdef"), tseq("gH") + dna_G,
                                         tseq("yz"))));
  }
}

INSTANTIATE_TEST_CASE_P(fwd_tracer_test, tracer_test,
                        ::testing::Values(std::make_pair(false /* not rev_comp */,
                                                         false /* tracer, not pop_tracer */)));
INSTANTIATE_TEST_CASE_P(rev_tracer_test, tracer_test,
                        ::testing::Values(std::make_pair(true /*rev_comp */,
                                                         false /* tracer, not pop_tracer */)));

INSTANTIATE_TEST_CASE_P(DISABLED_bidir_tracer_test, tracer_test,
                        ::testing::Values(std::make_pair(false /* not rev_comp */,
                                                         true /* pop_tracer */)));

}  // namespace variants

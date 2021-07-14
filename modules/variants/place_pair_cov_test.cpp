#include "modules/variants/place_pair_cov.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/pair_cov.h"
#include "modules/variants/place_pair_cov.h"
#include "modules/variants/read_cov.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

struct placed_read {
  // Sequence in starting assembly before the read starts
  dna_sequence before_read;

  // Sequence during the read
  dna_sequence read_seq;

  // Sequence in starting assembly after the read ends
  dna_sequence after_read;

  // Assembly ids the read traverses
  std::vector<int> assembly_ids;
};

dna_sequence repeat(dna_sequence seq, int n) {
  dna_sequence result;
  while (n) {
    result += seq;
    --n;
  }
  return result;
}

std::ostream& operator<<(std::ostream& os, const placed_read& pr) {
  return os << "\nPlaced(before=" << pr.before_read << ", read=" << pr.read_seq
            << ", after=" << pr.after_read << ", assembly_ids=" << PrintToString(pr.assembly_ids)
            << ")";
}

MATCHER_P4(_PlacedReadIs, before, seq, after, aids, "") {
  bool ok = true;
  if (before != arg.before_read) {
    (*result_listener) << " where the before seq " << arg.before_read << " is not the expected "
                       << before;
    ok = false;
  }
  if (seq != arg.read_seq) {
    (*result_listener) << " where the read seq " << arg.read_seq << " is not the expected " << seq;
    ok = false;
  }
  if (after != arg.after_read) {
    (*result_listener) << " where the after seq " << arg.after_read << " is not the expected "
                       << after;
    ok = false;
  }
  if (!ok) {
    (*result_listener) << "\n";
  }
  if (!ExplainMatchResult(ContainerEq(aids), arg.assembly_ids, result_listener)) {
    ok = false;
  }
  return ok;
}

Matcher<const placed_read&> PlacedReadIs(dna_sequence before, dna_sequence seq, dna_sequence after,
                                         std::vector<int> aids) {
  return _PlacedReadIs(before, seq, after, aids);
}

class place_pair_cov_test_base : public assemble_test {
 public:
  using placer_inspector_t = place_pair_cov::placer_inspector_t;

  void SetUp() override {
    m_options.min_pair_distance = 0;
    m_options.max_pair_distance = tseq("a").size() * 100;
  }
  void start_placer(aoffset_t ideal_len, placer_inspector_t inspect = placer_inspector_t()) {
    m_popts.ideal_pair_distance = ideal_len;
    auto place_step = make_unique<place_pair_cov>(m_options, m_popts, test_output());
    place_step->testing_set_inspector(inspect);
    auto read_step = make_unique<read_cov>(m_options, std::move(place_step));
    m_pipeline = std::move(read_step);
  }

  void add_ref_assembly(int asm_id, aoffset_t left_offset, aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = asm_id;

    a->tags.insert("place_pair_cov_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->matches_reference = true;
    a->seq = m_scaffold.subscaffold(left_offset, right_offset - left_offset).get_simple();

    m_pipeline->add(std::move(a));
  }

  void add_var_assembly(int asm_id, aoffset_t left_offset, dna_sequence seq,
                        aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = asm_id;

    a->tags.insert("place_pair_cov_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;

    m_pipeline->add(std::move(a));
  }

  void done() {
    m_pipeline->flush();
    m_pipeline.reset();

    constexpr bool k_dump_reads = false;
    read_id_set all_read_ids;
    for (const auto& a : m_assemblies) {
      all_read_ids |= a.pair_read_coverage->all_read_ids();
      if (k_dump_reads) {
        std::cerr << "Asm: " << a << "\n";
        for (const auto& cov_entry : a.pair_read_coverage->reads()) {
          std::cerr << "Pair cov entry: " << cov_entry
                    << " reads=" << PrintToString(cov_entry.read_ids.to_vector()) << "\n";
        }
        std::cerr << a.read_coverage->reads().size() << " read cov entries total\n";
        for (const auto& cov_entry : a.read_coverage->reads()) {
          std::cerr << "Read cov entry: " << cov_entry
                    << " reads=" << PrintToString(cov_entry.read_ids.to_vector()) << "\n";
        }
      }
    }

    for (uint32_t read_id : all_read_ids) {
      if (k_dump_reads) {
        std::cerr << "Tracking " << read_id << "\n";
      }
      placed_read pr;
      auto rd = m_options.readmap->get_read_by_id(read_id);
      int read_len = rd.size();
      pr.read_seq = rd.get_seqset_entry().sequence();
      dna_sequence read_seq_left = pr.read_seq;

      aoffset_t ref_offset = -1;
      aoffset_t read_offset = -1;

      for (const auto& a : m_assemblies) {
        if (k_dump_reads) {
          std::cerr << "Checking asm " << a << "\n";
        }
        for (const auto& cov_entry : a.pair_read_coverage->reads()) {
          if (k_dump_reads) {
            std::cerr << "Checking cov entry " << cov_entry << "\n";
          }
          if (!cov_entry.read_ids.contains(read_id)) {
            continue;
          }

          pr.assembly_ids.push_back(a.assembly_id);

          dna_sequence match_seq = a.seq;
          aoffset_t offset_from_end = -1;

          auto read_info = [&]() {
            std::stringstream out;

            out << "\n"
                << a << "\ncov_entry=" << cov_entry << " read_id=" << read_id
                << ": placed so far=" << pr << "\nread offset=" << read_offset
                << " offset from end=" << offset_from_end << " ref offset=" << ref_offset
                << " match seq=" << match_seq;
            return out.str();
          };

          if (cov_entry.offset >= 0) {
            // Starts in this assembly.
            EXPECT_EQ(ref_offset, -1);
            EXPECT_EQ(read_offset, -1);

            EXPECT_EQ(cov_entry.read_len, read_len) << read_info();
            EXPECT_EQ(pr.before_read, dna_sequence()) << read_info();
            EXPECT_EQ(pr.after_read, dna_sequence()) << read_info();

            pr.before_read = a.seq.subseq(0, cov_entry.offset);
            match_seq = match_seq.subseq(cov_entry.offset, match_seq.size() - cov_entry.offset);
            read_offset = 0;
          } else {
            EXPECT_EQ(ref_offset, a.left_offset) << read_info();
          }
          ref_offset = a.right_offset;

          offset_from_end = aoffset_t(a.seq.size()) - (cov_entry.offset + cov_entry.read_len);
          if (offset_from_end >= 0) {
            // Ends in this assembly
            EXPECT_LE(offset_from_end, match_seq.size()) << read_info();
            EXPECT_EQ(pr.after_read, dna_sequence()) << read_info();
            pr.after_read = a.seq.subseq(a.seq.size() - offset_from_end, offset_from_end);
            match_seq = match_seq.subseq(0, match_seq.size() - offset_from_end);

            EXPECT_EQ(read_offset + match_seq.size(), pr.read_seq.size()) << read_info();
          } else {
            EXPECT_LT(read_offset + match_seq.size(), pr.read_seq.size()) << read_info();
          }
          EXPECT_EQ(match_seq, pr.read_seq.subseq(read_offset, match_seq.size())) << read_info();
          read_offset += match_seq.size();
        }
      }

      if (rd.has_mate()) {
        EXPECT_EQ(read_len, read_offset)
            << ": placed so far=" << pr << " read offset=" << read_offset
            << " ref offset=" << ref_offset;
        m_placed_reads.emplace_back(std::move(pr));
      } else {
        EXPECT_EQ(read_offset, -1)
            << ": Should not have placed non-mated read!  placed so far=" << pr
            << " read offset=" << read_offset << " ref offset=" << ref_offset;
      }
    }
  }

  std::vector<placed_read> m_placed_reads;
  pipeline_step_t m_pipeline;
  bool m_skip_ambiguous = false;
  place_pair_options m_popts;
};

class place_pair_cov_test : public place_pair_cov_test_base, public WithParamInterface<bool> {
 public:
  void SetUp() override {
    m_skip_ambiguous = GetParam();
    if (m_skip_ambiguous) {
      m_popts.max_ambig = 1;
    }
    place_pair_cov_test_base::SetUp();
  }
};

TEST_P(place_pair_cov_test, simple_ref_single_assembly) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdef"), tseq_rc("vwxy")}});

  start_placer(100, [&](place_pair_cov* pl) {
    aoffset_t end_offset = tseq("abcdefghijklmnopqrstuvwxyz").size();
    pl->dump_state("Inspecting");
    EXPECT_THAT(pl->testing_distances_between(0, 0), ElementsAre(0));
    EXPECT_THAT(pl->testing_distances_between(0, end_offset), ElementsAre(end_offset));
    EXPECT_THAT(pl->testing_distances_between(end_offset, end_offset), ElementsAre(0));
  });

  add_ref_assembly(123, 0, tseq("abcdefghijklmnopqrstuvwxyz").size());
  done();

  EXPECT_THAT(
      m_placed_reads,
      ElementsAre(PlacedReadIs(tseq("a"), tseq("bcdef"), tseq("ghijklmnopqrstuvwxyz"), {123}),
                  PlacedReadIs(tseq("abcdefghijklmnopqrstu"), tseq("vwxy"), tseq("z"), {123})));
}

TEST_P(place_pair_cov_test, simple_var_multi_assembly_edges) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdeF"), tseq_rc("Vwxy")}});

  start_placer(121, [&](place_pair_cov* pl) {
    aoffset_t ref_a_offset = tseq("a").size();
    aoffset_t ref_y_offset = tseq("abcdefghijklmnopqrstuvwxy").size();
    aoffset_t var_start_offset = tseq("abcde").size();
    aoffset_t var_end_offset = tseq("abcdefghijklmnopqrstuv").size();

    aoffset_t ref_end_offset = tseq("abcdefghijklmnopqrstuvwxyz").size();

    EXPECT_THAT(pl->testing_distances_between(0, ref_a_offset), ElementsAre(tseq("a").size()));
    EXPECT_THAT(pl->testing_distances_between(ref_y_offset, ref_end_offset),
                ElementsAre(tseq("z").size()));

    EXPECT_THAT(pl->testing_distances_between(0, var_start_offset),
                ElementsAre(tseq("abcde").size()));
    EXPECT_THAT(pl->testing_distances_between(var_end_offset, ref_end_offset),
                ElementsAre(tseq("wxyz").size()));

    EXPECT_THAT(pl->testing_distances_between(var_start_offset, var_end_offset),
                ElementsAre(tseq("F_V").size(), tseq("fghijklmnopqrstuv").size()));
    EXPECT_THAT(pl->testing_distances_between(ref_a_offset, ref_y_offset),
                ElementsAre(tseq("bcdeF_Vwxy").size()
                            // This distance would be present if we
                            // kept more than one distance beyond
                            // ideal:
                            // tseq("bcdefghijklmnopqrstuvwxy").size()
                            ));
    EXPECT_THAT(pl->testing_distances_between(0, ref_end_offset),
                ElementsAre(tseq("abcdeF_Vwxyz").size()
                            // This distance would be present if we
                            // kept more than one distance beyond
                            // ideal:
                            //  tseq("abcdefghijklmnopqrstuvwxyz").size()
                            ));
  });

  add_ref_assembly(1, 0, tseq("a").size());
  add_ref_assembly(2, tseq("a").size(), tseq("abcde").size());
  add_ref_assembly(3, tseq("abcde").size(), tseq("abcdefghijklmnopqrstuv").size());
  add_var_assembly(4, tseq("abcde").size(), tseq("F_V"), tseq("abcdefghijklmnopqrstuv").size());
  add_ref_assembly(5, tseq("abcdefghijklmnopqrstuv").size(),
                   tseq("abcdefghijklmnopqrstuvwxy").size());
  add_ref_assembly(6, tseq("abcdefghijklmnopqrstuvwxy").size(),
                   tseq("abcdefghijklmnopqrstuvwxyz").size());

  done();

  EXPECT_THAT(m_placed_reads,
              ElementsAre(PlacedReadIs(tseq(""), tseq("bcdeF"), tseq("_V"), {2, 4}),
                          PlacedReadIs(tseq("F_"), tseq("Vwxy"), tseq(""), {4, 5})));
}

TEST_P(place_pair_cov_test, pick_best_choice_single_var) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdeF"), tseq_rc("Vwxy")}});

  aoffset_t ideal_dist = tseq("bcdeF_").size() + 1 + tseq("_Vwxy").size();

  start_placer(ideal_dist, [&](place_pair_cov* pl) {
    aoffset_t var_start_offset = tseq("abcde").size();
    aoffset_t var_end_offset = tseq("abcdefghijklmnopqrstuv").size();

    EXPECT_THAT(pl->testing_distances_between(var_start_offset, var_end_offset),
                ElementsAre(tseq("F__V").size(), tseq("F_").size() + 1 + tseq("_V").size(),
                            tseq("F_").size() + 2 + tseq("_V").size(),
                            // And the ref path is always available:
                            tseq("fghijklmnopqrstuv").size()));
  });

  add_ref_assembly(1, 0, tseq("abcde").size());
  // Shorter than ideal:
  add_var_assembly(2, tseq("abcde").size(), tseq("F__V"), tseq("abcdefghijklmnopqrstuv").size());
  // Just right:
  add_var_assembly(3, tseq("abcde").size(), tseq("F_") + dna_T + tseq("_V"),
                   tseq("abcdefghijklmnopqrstuv").size());
  // Longer than ideal:
  add_var_assembly(4, tseq("abcde").size(), tseq("F_") + dna_T + dna_T + tseq("_V"),
                   tseq("abcdefghijklmnopqrstuv").size());
  add_ref_assembly(5, tseq("abcdefghijklmnopqrstuv").size(),
                   tseq("abcdefghijklmnopqrstuvwxyz").size());

  done();

  EXPECT_THAT(
      m_placed_reads,
      ElementsAre(PlacedReadIs(tseq("a"), tseq("bcdeF"), tseq("_") + dna_T + tseq("_V"), {1, 3}),
                  PlacedReadIs(tseq("F_") + dna_T + tseq("_"), tseq("Vwxy"), tseq("z"), {3, 5})));
}

// Where the inside ends of the mates overlap, so like:
// Read1:      abcdefghijkl
// Read2(rc):          ijklmnopq
TEST_P(place_pair_cov_test, pick_best_choice_single_var_overlap) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads(
      {{tseq("bcdeF") + repeat(dna_A, 10), (repeat(dna_A, 10) + tseq("Mnopq")).rev_comp()}});

  aoffset_t ideal_dist = tseq("bcdeF").size() + 15 + tseq("Mnopq").size();

  start_placer(ideal_dist);

  add_ref_assembly(1, 0, tseq("abcde").size());
  // Shorter than ideal:
  add_var_assembly(2, tseq("abcde").size(), tseq("F") + repeat(dna_A, 14) + tseq("M"),
                   tseq("abcdefghijklm").size());
  // Just right:
  add_var_assembly(3, tseq("abcde").size(), tseq("F") + repeat(dna_A, 15) + tseq("M"),
                   tseq("abcdefghijklm").size());
  // Longer than ideal:
  add_var_assembly(4, tseq("abcde").size(), tseq("F") + repeat(dna_A, 16) + tseq("M"),
                   tseq("abcdefghijklm").size());

  add_ref_assembly(5, tseq("abcde").size(), tseq("abcdefghijklm").size());
  add_ref_assembly(6, tseq("abcdefghijklm").size(), tseq("abcdefghijklmnopqrstuvwxyz").size());

  done();

  EXPECT_THAT(m_placed_reads,
              ElementsAre(  //
                  PlacedReadIs(tseq("F") + repeat(dna_A, 5), repeat(dna_A, 10) + tseq("Mnopq"),
                               tseq("rstuvwxyz"), {3, 6}),
                  PlacedReadIs(tseq("a"), tseq("bcdeF") + repeat(dna_A, 10),
                               repeat(dna_A, 5) + tseq("M"), {1, 3})));
}

TEST_P(place_pair_cov_test, pick_best_choice_multi_var_overlap) {
  use_ref_parts({{0, tseq("abcdef") + repeat(dna_A, 10) + tseq("ghijkl")}});
  aoffset_t start_g = tseq("abcdef").size() + 10;
  use_paired_reads(
      {{tseq("bcdeF") + repeat(dna_A, 15), (repeat(dna_A, 15) + tseq("ghijk")).rev_comp()}});

  aoffset_t ideal_dist = tseq("bcdeF").size() + 20 + tseq("ghijk").size();

  start_placer(ideal_dist);

  add_ref_assembly(1, 0, tseq("abcde").size());
  // Shorter than ideal:
  add_var_assembly(2, tseq("abcde").size(), tseq("F") + repeat(dna_A, 9), tseq("abcdef").size());
  // Just right:
  add_var_assembly(3, tseq("abcde").size(), tseq("F") + repeat(dna_A, 10), tseq("abcdef").size());
  // Longer than ideal:
  add_var_assembly(4, tseq("abcde").size(), tseq("F") + repeat(dna_A, 11), tseq("abcdef").size());

  add_ref_assembly(5, tseq("abcde").size(), tseq("abcdef").size());
  add_ref_assembly(6, tseq("abcdef").size(), start_g + tseq("ghijkl").size());

  done();

  EXPECT_THAT(m_placed_reads,
              ElementsAre(  //
                  PlacedReadIs(tseq("F") + repeat(dna_A, 5), repeat(dna_A, 15) + tseq("ghijk"),
                               tseq("l"), {3, 6}),
                  PlacedReadIs(tseq("a"), tseq("bcdeF") + repeat(dna_A, 15),
                               repeat(dna_A, 5) + tseq("ghijkl"), {1, 3, 6})));
}

TEST_P(place_pair_cov_test, pick_best_choice_multi_insert) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdef"), tseq_rc("MATE")}});

  aoffset_t ideal_dist = tseq("bcdefghijk").size() + 1 + tseq("MATE").size();

  start_placer(ideal_dist);

  add_ref_assembly(1, 0, tseq("abcdef").size());
  // Shorter than ideal:
  add_var_assembly(2, tseq("abcdefghijk").size(), tseq("MATE"), tseq("abcdefghijk").size());
  // Just right:
  add_var_assembly(3, tseq("abcdefghijk").size() + 1, tseq("MATE"), tseq("abcdefghijk").size() + 1);
  // Longer than ideal:
  add_var_assembly(4, tseq("abcdefghijk").size() + 2, tseq("MATE"), tseq("abcdefghijk").size() + 2);

  done();

  EXPECT_THAT(m_placed_reads, ElementsAre(PlacedReadIs(tseq("a"), tseq("bcdef"), tseq(""), {1}),
                                          PlacedReadIs(tseq(""), tseq("MATE"), tseq(""), {3})));
}

TEST_P(place_pair_cov_test, round_robin_between_equals) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("bcdef"), tseq_rc("MATE")},
                    {tseq("bcdef"), tseq_rc("MATE")},
                    {tseq("bcdef"), tseq_rc("MATE")},
                    {tseq("bcdef"), tseq_rc("MATE")}});

  aoffset_t ideal_dist = tseq("bcdefghijk").size() + 1 + tseq("MATE").size();

  start_placer(ideal_dist);

  add_ref_assembly(1, 0, tseq("abcdef").size());
  add_var_assembly(2, tseq("abcdefghijk").size() + 1, tseq("MATE"), tseq("abcdefghijk").size() + 1);
  add_var_assembly(3, tseq("abcdefghijk").size(), dna_G + tseq("MATE"), tseq("abcdefghijk").size());

  done();

  if (m_skip_ambiguous) {
    EXPECT_THAT(m_placed_reads, IsEmpty());
  } else {
    EXPECT_THAT(m_placed_reads,
                UnorderedElementsAre(  //

                    // forward reads:
                    PlacedReadIs(tseq("a"), tseq("bcdef"), tseq(""), {1}),
                    PlacedReadIs(tseq("a"), tseq("bcdef"), tseq(""), {1}),
                    PlacedReadIs(tseq("a"), tseq("bcdef"), tseq(""), {1}),
                    PlacedReadIs(tseq("a"), tseq("bcdef"), tseq(""), {1}),

                    // mate reads:
                    PlacedReadIs(tseq(""), tseq("MATE"), tseq(""), {2}),
                    PlacedReadIs(tseq(""), tseq("MATE"), tseq(""), {2}),

                    PlacedReadIs(dna_G, tseq("MATE"), tseq(""), {3}),
                    PlacedReadIs(dna_G, tseq("MATE"), tseq(""), {3})));
  }
}

TEST_P(place_pair_cov_test, round_robin_trace_out) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads({{tseq("READf"), tseq_rc("kMATE")},
                    {tseq("READf"), tseq_rc("kMATE")},
                    {tseq("READf"), tseq_rc("kMATE")},
                    {tseq("READf"), tseq_rc("kMATE")}});

  aoffset_t ideal_dist = tseq("bcdefghijkMATE").size();

  start_placer(ideal_dist);

  add_var_assembly(1, 0, tseq("READ"), tseq("abcde").size());
  add_var_assembly(2, 0, tseq("READ"), tseq("abcde").size());

  add_ref_assembly(3, tseq("abcde").size(), tseq("abcdefghijk").size());

  add_var_assembly(4, tseq("abcdefghijk").size(), tseq("MATE"), 150);
  add_var_assembly(5, tseq("abcdefghijk").size(), tseq("MATE"), 160);

  done();

  EXPECT_THAT(m_placed_reads,
              UnorderedElementsAre(  //
                                     // forward reads:
                  PlacedReadIs(tseq(""), tseq("READf"), tseq("ghijk"), {1, 3}),
                  PlacedReadIs(tseq(""), tseq("READf"), tseq("ghijk"), {1, 3}),
                  PlacedReadIs(tseq(""), tseq("READf"), tseq("ghijk"), {2, 3}),
                  PlacedReadIs(tseq(""), tseq("READf"), tseq("ghijk"), {2, 3}),

                  // mate reads:
                  PlacedReadIs(tseq("fghij"), tseq("kMATE"), tseq(""), {3, 4}),
                  PlacedReadIs(tseq("fghij"), tseq("kMATE"), tseq(""), {3, 4}),
                  PlacedReadIs(tseq("fghij"), tseq("kMATE"), tseq(""), {3, 5}),
                  PlacedReadIs(tseq("fghij"), tseq("kMATE"), tseq(""), {3, 5})));
}

INSTANTIATE_TEST_CASE_P(place_pair_cov_tests, place_pair_cov_test,
                        ::testing::Values(true, false) /* skip ambiguous */);

class place_pair_cov_bounds_test
    : public place_pair_cov_test_base,
      public testing::WithParamInterface<
          std::pair<bool /* decrement minimum */, bool /* increment maximum */>> {
 public:
  void SetUp() override {
    // m_dec_min and m_inc_max are used to make sure tests fail if
    // they're not correctly testing bounds enforcement.
    std::tie(m_dec_min, m_inc_max) = GetParam();
    place_pair_cov_test_base::SetUp();
  }

  void start_placer(aoffset_t ideal_len, placer_inspector_t inspect = placer_inspector_t()) {
    if (m_dec_min) {
      --m_options.min_pair_distance;
    } else if (m_inc_max) {
      ++m_options.max_pair_distance;
    }
    place_pair_cov_test_base::start_placer(ideal_len, inspect);
  }

  Matcher<const std::vector<placed_read>&> BoundsCheck(
      const Matcher<const std::vector<placed_read>&>& m) {
    if (m_dec_min || m_inc_max) {
      return Not(m);
    } else {
      return m;
    }
  }

  bool m_dec_min, m_inc_max;
};

TEST_P(place_pair_cov_bounds_test, single_asm_bounds) {
  use_ref_parts({{0, tseq("abcde") + dna_A + dna_T + tseq("fghijkl") + dna_A + dna_T +
                         tseq("mnopqrstuvwxyz")}});

  use_paired_reads({// Correct distance:
                    {dna_T + tseq("fgh"), (tseq("jkl") + dna_A).rev_comp()},
                    // Too big:
                    {dna_A + dna_T + tseq("fgh"), (tseq("jkl") + dna_A).rev_comp()},
                    {dna_T + tseq("fgh"), (tseq("jkl") + dna_A + dna_T).rev_comp()},
                    // Too small:
                    {tseq("fgh"), (tseq("jkl") + dna_A).rev_comp()},
                    {dna_T + tseq("fgh"), (tseq("jkl")).rev_comp()}});
  aoffset_t valid_dist = (dna_T + tseq("fghijkl") + dna_A).size();
  m_options.min_pair_distance = valid_dist;
  m_options.max_pair_distance = valid_dist;
  start_placer(valid_dist);
  add_ref_assembly(
      1, 0, tseq("abcde").size() + 2 + tseq("fghijkl").size() + 2 + tseq("mnopqrstuvwxyz").size());

  done();

  EXPECT_THAT(m_placed_reads,
              BoundsCheck(UnorderedElementsAre(  //
                  PlacedReadIs(tseq("abcde") + dna_A, dna_T + tseq("fgh"),
                               tseq("ijkl") + dna_A + dna_T + tseq("mnopqrstuvwxyz"), {1}),
                  PlacedReadIs(tseq("abcde") + dna_A + dna_T + tseq("fghi"), tseq("jkl") + dna_A,
                               dna_T + tseq("mnopqrstuvwxyz"), {1}))));
}

TEST_P(place_pair_cov_bounds_test, single_asm_bounds_overlap) {
  use_ref_parts({{0, tseq("abcde") + dna_A + dna_T + tseq("fghijkl") + dna_A + dna_T +
                         tseq("mnopqrstuvwxyz")}});

  use_paired_reads({// Correct distance:
                    {dna_T + tseq("fghi"), (tseq("ijkl") + dna_A).rev_comp()},
                    // Too big:
                    {dna_A + dna_T + tseq("fghi"), (tseq("ijkl") + dna_A).rev_comp()},
                    {dna_T + tseq("fghi"), (tseq("ijkl") + dna_A + dna_T).rev_comp()},
                    // Too small:
                    {tseq("fghi"), (tseq("ijkl") + dna_A).rev_comp()},
                    {dna_T + tseq("fghi"), (tseq("ijkl")).rev_comp()}});
  aoffset_t valid_dist = (dna_T + tseq("fghijkl") + dna_A).size();
  m_options.min_pair_distance = valid_dist;
  m_options.max_pair_distance = valid_dist;
  start_placer(valid_dist);
  add_ref_assembly(
      1, 0, tseq("abcde").size() + 2 + tseq("fghijkl").size() + 2 + tseq("mnopqrstuvwxyz").size());

  done();

  EXPECT_THAT(m_placed_reads,
              BoundsCheck(UnorderedElementsAre(  //
                  PlacedReadIs(tseq("abcde") + dna_A, dna_T + tseq("fghi"),
                               tseq("jkl") + dna_A + dna_T + tseq("mnopqrstuvwxyz"), {1}),
                  PlacedReadIs(tseq("abcde") + dna_A + dna_T + tseq("fgh"), tseq("ijkl") + dna_A,
                               dna_T + tseq("mnopqrstuvwxyz"), {1}))));
}

TEST_P(place_pair_cov_bounds_test, two_asm_bounds) {
  use_ref_parts({{0, tseq("abcde") + dna_A + dna_T + tseq("fghijklm") + dna_A + dna_T +
                         tseq("nopqrstuvwxyz")}});

  use_paired_reads({// Correct distance:
                    {dna_T + tseq("fghi"), (tseq("jklm") + dna_A).rev_comp()},
                    // Too big:
                    {dna_A + dna_T + tseq("fghi"), (tseq("jklm") + dna_A).rev_comp()},
                    {dna_T + tseq("fghi"), (tseq("jklm") + dna_A + dna_T).rev_comp()},
                    // Too small:
                    {tseq("fghi"), (tseq("jklm") + dna_A).rev_comp()},
                    {dna_T + tseq("fghi"), (tseq("jklm")).rev_comp()}});
  aoffset_t valid_dist = (dna_T + tseq("fghijklm") + dna_A).size();
  m_options.min_pair_distance = valid_dist;
  m_options.max_pair_distance = valid_dist;
  start_placer(valid_dist);
  aoffset_t asm_split_point = tseq("abcde").size() + 2 + tseq("fghi").size();
  add_ref_assembly(1, 0, asm_split_point);
  add_ref_assembly(2, asm_split_point,
                   asm_split_point + tseq("jklm").size() + 2 + tseq("nopqrstuvwxyz").size());

  done();

  EXPECT_THAT(m_placed_reads,
              BoundsCheck(UnorderedElementsAre(  //
                  PlacedReadIs(tseq("abcde") + dna_A, dna_T + tseq("fghi"), {}, {1}),
                  PlacedReadIs({}, tseq("jklm") + dna_A, dna_T + tseq("nopqrstuvwxyz"), {2}))));
}

TEST_P(place_pair_cov_bounds_test, two_asm_bounds_overlap) {
  use_ref_parts({{0, tseq("abcde") + dna_A + dna_T + tseq("fghijklm") + dna_A + dna_T +
                         tseq("nopqrstuvwxyz")}});

  use_paired_reads({// Correct distance:
                    {dna_T + tseq("fghij"), (tseq("ijklm") + dna_A).rev_comp()},
                    // Too big:
                    {dna_A + dna_T + tseq("fghij"), (tseq("ijklm") + dna_A).rev_comp()},
                    {dna_T + tseq("fghij"), (tseq("ijklm") + dna_A + dna_T).rev_comp()},
                    // Too small:
                    {tseq("fghij"), (tseq("ijklm") + dna_A).rev_comp()},
                    {dna_T + tseq("fghij"), (tseq("ijklm")).rev_comp()}});
  aoffset_t valid_dist = (dna_T + tseq("fghijklm") + dna_A).size();
  m_options.min_pair_distance = valid_dist;
  m_options.max_pair_distance = valid_dist;
  start_placer(valid_dist);
  aoffset_t asm_split_point = tseq("abcde").size() + 2 + tseq("fghi").size();
  add_ref_assembly(1, 0, asm_split_point);
  add_ref_assembly(2, asm_split_point,
                   asm_split_point + tseq("jklm").size() + 2 + tseq("nopqrstuvwxyz").size());

  done();

  EXPECT_THAT(m_placed_reads,
              BoundsCheck(UnorderedElementsAre(  //
                  PlacedReadIs(tseq("abcde") + dna_A, dna_T + tseq("fghij"),
                               tseq("klm") + dna_A + dna_T + tseq("nopqrstuvwxyz"), {1, 2}),
                  PlacedReadIs(tseq("abcde") + dna_A + dna_T + tseq("fgh"), tseq("ijklm") + dna_A,
                               dna_T + tseq("nopqrstuvwxyz"), {1, 2}))));
}

TEST_P(place_pair_cov_bounds_test, three_asm_bounds) {
  use_ref_parts({{0, tseq("abcde") + dna_A + dna_T + tseq("fghi") + dna_G + tseq("jklm") + dna_A +
                         dna_T + tseq("nopqrstuvwxyz")}});

  use_paired_reads({// Correct distance:
                    {dna_T + tseq("fgh"), (tseq("klm") + dna_A).rev_comp()},
                    // Too big:
                    {dna_A + dna_T + tseq("fgh"), (tseq("klm") + dna_A).rev_comp()},
                    {dna_T + tseq("fgh"), (tseq("klm") + dna_A + dna_T).rev_comp()},
                    // Too small:
                    {tseq("fgh"), (tseq("klm") + dna_A).rev_comp()},
                    {dna_T + tseq("fgh"), (tseq("klm")).rev_comp()}});
  aoffset_t valid_dist = (dna_T + tseq("fghi") + dna_G + tseq("jklm") + dna_A).size();
  m_options.min_pair_distance = valid_dist;
  m_options.max_pair_distance = valid_dist;
  start_placer(valid_dist);
  aoffset_t asm_split_point1 = tseq("abcde").size() + 2 + tseq("fghi").size();
  add_ref_assembly(1, 0, asm_split_point1);
  aoffset_t asm_split_point2 = asm_split_point1 + 1;
  add_ref_assembly(2, asm_split_point1, asm_split_point2);
  add_ref_assembly(3, asm_split_point2,
                   asm_split_point2 + tseq("jklm").size() + 2 + tseq("nopqrstuvwxyz").size());

  done();

  EXPECT_THAT(
      m_placed_reads,
      BoundsCheck(UnorderedElementsAre(  //
          PlacedReadIs(tseq("abcde") + dna_A, dna_T + tseq("fgh"), tseq("i"), {1}),
          PlacedReadIs(tseq("j"), tseq("klm") + dna_A, dna_T + tseq("nopqrstuvwxyz"), {3}))));
}

TEST_P(place_pair_cov_bounds_test, three_asm_bounds_overlap) {
  use_ref_parts({{0, tseq("abcde") + dna_A + dna_T + tseq("fghi") + dna_G + tseq("jklm") + dna_A +
                         dna_T + tseq("nopqrstuvwxyz")}});

  use_paired_reads(
      {// Correct distance:
       {dna_T + tseq("fghi") + dna_G + tseq("j"),
        (tseq("i") + dna_G + tseq("jklm") + dna_A).rev_comp()},
       // Too big:
       {dna_A + dna_T + tseq("fghi") + dna_G + tseq("j"),
        (tseq("i") + dna_G + tseq("jklm") + dna_A).rev_comp()},
       {dna_T + tseq("fghi") + dna_G + tseq("j"),
        (tseq("i") + dna_G + tseq("jklm") + dna_A + dna_T).rev_comp()},
       // Too small:
       {tseq("fghi") + dna_G + tseq("j"), (tseq("i") + dna_G + tseq("jklm") + dna_A).rev_comp()},
       {dna_T + tseq("fghi") + dna_G + tseq("j"), (tseq("i") + dna_G + tseq("jklm")).rev_comp()}});
  aoffset_t valid_dist = (dna_T + tseq("fghi") + dna_G + tseq("jklm") + dna_A).size();
  m_options.min_pair_distance = valid_dist;
  m_options.max_pair_distance = valid_dist;
  start_placer(valid_dist);
  aoffset_t asm_split_point1 = tseq("abcde").size() + 2 + tseq("fghi").size();
  add_ref_assembly(1, 0, asm_split_point1);
  aoffset_t asm_split_point2 = asm_split_point1 + 1;
  add_ref_assembly(2, asm_split_point1, asm_split_point2);
  add_ref_assembly(3, asm_split_point2,
                   asm_split_point2 + tseq("jklm").size() + 2 + tseq("nopqrstuvwxyz").size());

  done();

  EXPECT_THAT(m_placed_reads,
              BoundsCheck(UnorderedElementsAre(  //
                  PlacedReadIs(tseq("abcde") + dna_A, dna_T + tseq("fghi") + dna_G + tseq("j"),
                               tseq("klm") + dna_A + dna_T + tseq("nopqrstuvwxyz"), {1, 2, 3}),
                  PlacedReadIs(tseq("abcde") + dna_A + dna_T + tseq("fgh"),
                               tseq("i") + dna_G + tseq("jklm") + dna_A,
                               dna_T + tseq("nopqrstuvwxyz"), {1, 2, 3}))));
}

INSTANTIATE_TEST_CASE_P(place_pair_cov_bounds_tests, place_pair_cov_bounds_test,
                        ::testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                          std::make_pair(false, true)));

}  // namespace variants

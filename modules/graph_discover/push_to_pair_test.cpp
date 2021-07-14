#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/graph_discover/push_to_pair.h"
#include "modules/graph_discover/update_rc_seqset_entries.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/ref_map.h"

using namespace testing;
using namespace dna_testutil;

namespace variants {

class push_to_pair_discover_test : public assemble_test {
 public:
  const std::string k_tag = "push_to_pair";

  void start() {  //
    auto update = make_unique<update_rc_seqset_entries>(m_options, test_output());
    update->enable_self_test();
    auto discover = make_unique<push_to_pair_discover>(m_options, k_tag, std::move(update));
    m_discover = make_unique<update_rc_seqset_entries>(m_options, std::move(discover));
    m_discover->enable_self_test();
  }

  void flush() {
    m_discover->flush();
    EXPECT_TRUE(m_discover->self_test_succeeded());
    m_discover.reset();
  }

  void add_ref_asm(aoffset_t left_offset, dna_sequence seq) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->tags.insert("push_to_pair_discover_test");
    a->left_offset = left_offset;
    a->right_offset = left_offset + seq.size();
    a->seq = seq;
    a->matches_reference = true;

    m_discover->add(std::move(a));
  }

  void add_var_asm(optional_aoffset left_offset, dna_sequence seq, optional_aoffset right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->tags.insert("push_to_pair_discover_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    a->matches_reference = false;

    m_discover->add(std::move(a));
  }

  std::unique_ptr<update_rc_seqset_entries> m_discover;
  std::vector<assembly> m_traced;
};

TEST_F(push_to_pair_discover_test, simple_ref_only_no_branches) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcdefg"), tseq("defghij"), tseq("efghijklm")});

  start();
  add_ref_asm(0, tseq("abcdefghijklm"));
  flush();

  EXPECT_THAT(m_tag_assemblies[k_tag], IsEmpty());
}

TEST_F(push_to_pair_discover_test, ref_rejoin) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcdefg"), tseq("defghij"), tseq("efghijklm"), tseq("ef") + dna_G + tseq("gh"),
             tseq("ghijklmno")});

  m_options.min_overlap = tseq("ef").size();
  start();
  add_var_asm(tseq("abcd").size(), tseq("ef") + dna_G, optional_aoffset::none);
  add_ref_asm(tseq("abcde").size(), tseq("fghi"));
  flush();

  // May discover it multiple times as it continues tracing to get a
  // full anchor.
  auto anchor_match =
      AssemblyIs(tseq("abcd").size(), tseq("ef") + dna_G + tseq("ghi"), tseq("abcdefghi").size());
  auto continue_trace_match =
      AssemblyIs(tseq("abcd").size(), tseq("ef") + dna_G + tseq("ghi"), optional_aoffset::none);
  EXPECT_THAT(m_tag_assemblies[k_tag], Contains(anchor_match));
  EXPECT_THAT(m_tag_assemblies[k_tag], Contains(continue_trace_match));
  EXPECT_THAT(m_tag_assemblies[k_tag], Each(AnyOf(anchor_match, continue_trace_match)));
}

TEST_F(push_to_pair_discover_test, pair_extend_prev) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads(
      // Pair support:
      {{tseq("abc"), (dna_G + tseq("ghi")).rev_comp()}},
      // Unpaired path:
      {tseq("cdef") + dna_G + tseq("gh")});
  m_options.min_overlap = tseq("ef").size();

  start();
  add_ref_asm(0, tseq("abcd"));
  add_var_asm(tseq("abcd").size(), tseq("ef") + dna_G, optional_aoffset::none);
  flush();

  EXPECT_THAT(m_tag_assemblies[k_tag],
              ElementsAre(AssemblyIs(tseq("abcd").size(), tseq("ef") + dna_G + tseq("ghi"),
                                     optional_aoffset::none)));
}

TEST_F(push_to_pair_discover_test, pair_extend_other_assembly) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads(
      // Pair support:
      {{tseq("abcde"), (dna_G + tseq("ghi")).rev_comp()}},
      // Unpaired path:
      {tseq("cdef") + dna_G + tseq("gh")});
  m_options.min_overlap = tseq("ef").size();

  start();
  add_ref_asm(0, tseq("abcd"));
  add_var_asm(tseq("abcd").size(), tseq("ef") + dna_G, optional_aoffset::none);
  flush();

  EXPECT_THAT(m_tag_assemblies[k_tag],
              ElementsAre(AssemblyIs(tseq("abcd").size(), tseq("ef") + dna_G + tseq("ghi"),
                                     optional_aoffset::none)));
}

TEST_F(push_to_pair_discover_test, pair_extend_within_max_pair_distance) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads(
      // Pair support:
      {{tseq("abc"), (dna_G + tseq("uvw")).rev_comp()}},
      // Unpaired path:
      {tseq("pqrst") + dna_G + tseq("uv")});
  m_options.min_overlap = tseq("ef").size();
  m_options.max_pair_distance = tseq("abcdefghijklmnopqrst").size();

  start();
  add_ref_asm(0, tseq("abcd"));
  add_ref_asm(tseq("abcd").size(), tseq("efghijklmnopqr"));
  add_var_asm(tseq("abcdefghijklmnopqr").size(), tseq("st") + dna_G, optional_aoffset::none);
  flush();

  EXPECT_THAT(m_tag_assemblies[k_tag],
              ElementsAre(AssemblyIs(tseq("abcdefghijklmnopqr").size(),
                                     tseq("st") + dna_G + tseq("uvw"), optional_aoffset::none)));
}

TEST_F(push_to_pair_discover_test, DISABLED_pair_extend_outside_max_pair_distance) {
  // Make sure pair support goes away eventually.
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_paired_reads(
      // Pair support:
      {{tseq("abc"), (dna_G + tseq("uvw")).rev_comp()}},
      // Unpaired path:
      {tseq("pqrst") + dna_G + tseq("uv")});
  m_options.min_overlap = tseq("ef").size();
  m_options.max_pair_distance = tseq("abcdefghijkl").size();
  start();
  add_ref_asm(0, tseq("abcd"));
  add_ref_asm(tseq("abcd").size(), tseq("efghijklmnopqr"));
  add_var_asm(tseq("abcdefghijklmnopqr").size(), tseq("st") + dna_G, optional_aoffset::none);
  flush();

  EXPECT_THAT(m_tag_assemblies[k_tag], IsEmpty());
}

}  // namespace variants

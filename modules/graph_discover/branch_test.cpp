#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/graph_discover/branch.h"
#include "modules/graph_discover/update_rc_seqset_entries.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/ref_map.h"

using namespace testing;
using namespace dna_testutil;

namespace variants {

class branch_discover_test : public assemble_test {
 public:
  const std::string k_tag = "branch_discover";

  void start() {  //
    auto update = make_unique<update_rc_seqset_entries>(m_options, test_output());
    update->enable_self_test();
    auto discover = make_unique<branch_discover>(m_options, k_tag, std::move(update));
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
    a->tags.insert("branch_discover_test");
    a->left_offset = left_offset;
    a->right_offset = left_offset + seq.size();
    a->seq = seq;
    a->matches_reference = true;

    m_discover->add(std::move(a));
  }

  void add_var_asm(optional_aoffset left_offset, dna_sequence seq, optional_aoffset right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->tags.insert("branch_discover_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    a->matches_reference = false;

    m_discover->add(std::move(a));
  }

  std::unique_ptr<update_rc_seqset_entries> m_discover;
  std::vector<assembly> m_traced;
};

TEST_F(branch_discover_test, simple_ref_only_no_branches) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcdefg"), tseq("defghij"), tseq("efghijklm")});

  start();
  add_ref_asm(0, tseq("abcdefghijklm"));
  flush();

  EXPECT_THAT(m_tag_assemblies[k_tag], IsEmpty());
}

TEST_F(branch_discover_test, simple_ref_only_single_branch) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcdefg"), tseq("defghij"), tseq("efghijklm"), tseq("ef") + dna_G + tseq("gh")});

  m_options.min_overlap = tseq("ef").size();
  start();
  add_ref_asm(tseq("a").size(), tseq("bcdefghijklm"));
  flush();

  EXPECT_THAT(
      m_tag_assemblies[k_tag],
      ElementsAre(AssemblyIs(tseq("a").size(), tseq("bcdef") + dna_G, optional_aoffset::none)));
}

TEST_F(branch_discover_test, overlap_too_small) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcdefg"), tseq("defghij"), tseq("efghijklm"), tseq("ef") + dna_G + tseq("gh")});

  m_options.min_overlap = tseq("ef").size() + 1;
  start();
  add_ref_asm(tseq("a").size(), tseq("bcdefghijklm"));
  flush();

  EXPECT_THAT(m_tag_assemblies[k_tag], IsEmpty());
}

}  // namespace variants

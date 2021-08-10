#include "modules/variants/align_count.h"

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

MATCHER_P4(AlignCountIs, aid, read_lens, local_aligned, tot_aligned, "") {
  bool ok = true;
  if (size_t(aid) != arg.assembly_id) {
    (*result_listener) << " where the assembly id " << arg.assembly_id << " is not the expected "
                       << aid;
    ok = false;
  }
  if (size_t(read_lens) != arg.align_count->local_read_lens) {
    (*result_listener) << " where the read lens " << arg.align_count->local_read_lens
                       << " is not the expected " << read_lens;
    ok = false;
  }
  if (size_t(local_aligned) != arg.align_count->local_aligned_bases) {
    (*result_listener) << " where the total aligned bases " << arg.align_count->local_aligned_bases
                       << " is not the expected " << local_aligned;
    ok = false;
  }
  if (size_t(tot_aligned) != arg.align_count->tot_aligned_bases) {
    (*result_listener) << " where the total aligned bases " << arg.align_count->tot_aligned_bases
                       << " is not the expected " << tot_aligned;
    ok = false;
  }
  if (!ok) {
    (*result_listener) << "\n";
  }
  return ok;
}

class align_count_test : public assemble_test {
 public:
  void start_counter() {
    auto count_step = make_unique<align_count>(m_options, test_output());
    auto read_step = make_unique<read_cov>(m_options, std::move(count_step));
    m_pipeline = std::move(read_step);
  }

  void add_ref_assembly(int asm_id, aoffset_t left_offset, aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = asm_id;

    a->tags.insert("align_count_test");
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

    a->tags.insert("align_count_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;

    m_pipeline->add(std::move(a));
  }

  void done() {
    m_pipeline->flush();
    m_pipeline.reset();
  }

  pipeline_step_t m_pipeline;
};

TEST_F(align_count_test, simple_ref_single_assembly) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcde")});
  start_counter();
  add_ref_assembly(123, 0, tseq("abcdefghijklmnopqrstuvwxyz").size());
  done();

  EXPECT_THAT(m_assemblies, ElementsAre(AlignCountIs(123, tseq("abcde").size(),
                                                     tseq("abcde").size(), tseq("abcde").size())));
}

TEST_F(align_count_test, simple_ref_two_assembly) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcde")});
  start_counter();
  add_var_assembly(1, 0, tseq("abcde"), tseq("abcde").size());
  add_ref_assembly(2, 0, tseq("abcde").size());
  add_ref_assembly(3, tseq("abcde").size(), tseq("abcdefghi").size());
  done();

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(
          AlignCountIs(1, tseq("abcde").size(), tseq("abcde").size(), tseq("abcde").size() * 2),
          AlignCountIs(2, tseq("abcde").size(), tseq("abcde").size(), tseq("abcde").size() * 2),
          AlignCountIs(3, 0, 0, 0)));
}

TEST_F(align_count_test, simple_ref_two_assembly_partial_read) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcdefgh")});
  start_counter();
  add_var_assembly(1, 0, tseq("abcde"), tseq("abcde").size());
  add_ref_assembly(2, 0, tseq("abcde").size());
  add_ref_assembly(3, tseq("abcde").size(), tseq("abcdefghi").size());
  done();

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(
          AlignCountIs(1, tseq("abcdefgh").size(), tseq("abcde").size(), tseq("abcde").size() * 2),
          AlignCountIs(2, tseq("abcdefgh").size(), tseq("abcde").size(), tseq("abcde").size() * 2),
          AlignCountIs(3, tseq("abcdefgh").size(), tseq("fgh").size(), tseq("fgh").size())));
}

TEST_F(align_count_test, overlapping) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("ABCDE")});

  start_counter();
  add_var_assembly(1, 0, tseq("ABCDE"), 20);
  add_ref_assembly(2, 0, 10);
  add_var_assembly(3, 10, tseq("ABCDE"), 30);
  add_ref_assembly(4, 10, 20);
  add_var_assembly(5, 20, tseq("ABCDE"), 40);
  add_ref_assembly(6, 20, 30);
  add_ref_assembly(7, 30, 40);
  done();

  EXPECT_THAT(
      m_assemblies,
      ElementsAre(                                                                                //
          AlignCountIs(1, tseq("ABCDE").size(), tseq("ABCDE").size(), tseq("ABCDE").size() * 2),  //
          AlignCountIs(2, 0, 0, 0),                                                               //
          AlignCountIs(3, tseq("ABCDE").size(), tseq("ABCDE").size(), tseq("ABCDE").size() * 3),  //
          AlignCountIs(4, 0, 0, 0),                                                               //
          AlignCountIs(5, tseq("ABCDE").size(), tseq("ABCDE").size(), tseq("ABCDE").size() * 2),  //
          AlignCountIs(6, 0, 0, 0),                                                               //
          AlignCountIs(7, 0, 0, 0)));
}

}  // namespace variants

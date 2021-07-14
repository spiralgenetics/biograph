#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/graph_discover/discover.h"
#include "modules/graph_discover/update_rc_seqset_entries.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/ref_map.h"

using namespace testing;
using namespace dna_testutil;

namespace variants {

class discover_test : public assemble_test {
 public:
  class stub_discover : public graph_discover {
   public:
    stub_discover(discover_test* fixture, const assemble_options& options, pipeline_step_t output)
        : graph_discover(options, std::move(output)), m_fixture(fixture) {}

    void on_trace(const active_assembly* act) override {  //
      stub_callback_t cb = boost::any_cast<stub_callback_t>(act->a->user_data);
      cb(act);
      m_fixture->m_traced.push_back(*act->a);
    }

    discover_test* m_fixture = nullptr;
  };
  using active_assembly = graph_discover::active_assembly;
  using stub_callback_t = std::function<void(const active_assembly*)>;
  static void no_stub_callback(const active_assembly*) {}

  void start() {  //
    auto update = make_unique<update_rc_seqset_entries>(m_options, test_output());
    update->enable_self_test();
    auto discover = make_unique<stub_discover>(this, m_options, std::move(update));
    m_discover = make_unique<update_rc_seqset_entries>(m_options, std::move(discover));
    m_discover->enable_self_test();
  }

  void flush() {
    m_discover->flush();
    EXPECT_TRUE(m_discover->self_test_succeeded());
    m_discover.reset();
  }

  void add_ref_asm(aoffset_t left_offset, dna_sequence seq, stub_callback_t cb = no_stub_callback) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->tags.insert("discover_test");
    a->left_offset = left_offset;
    a->right_offset = left_offset + seq.size();
    a->seq = seq;
    a->matches_reference = true;
    a->user_data = cb;

    m_discover->add(std::move(a));
  }

  void add_var_asm(optional_aoffset left_offset, dna_sequence seq, optional_aoffset right_offset,
                   stub_callback_t cb = no_stub_callback) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->tags.insert("discover_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    a->matches_reference = false;
    a->user_data = cb;

    m_discover->add(std::move(a));
  }

  std::unique_ptr<update_rc_seqset_entries> m_discover;
  std::vector<assembly> m_traced;
};

TEST_F(discover_test, simple_ref_only) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcdefg"), tseq("defghij"), tseq("efghijklm")});

  start();
  add_ref_asm(0, tseq("abcdefghijklm"), [](const active_assembly* act) {
    EXPECT_THAT(act->a->rc_seqset_entries.ends(), ElementsAre(SeqsetEntryIs(dna_sequence())));
    EXPECT_THAT(act->a->rc_seqset_entries.starts(),
                ElementsAre(SeqsetEntryIs(tseq_rc("efghijklm"))));
  });
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(RefAssemblyIs(0, tseq("abcdefghijklm").size())));
  EXPECT_THAT(m_traced, ElementsAre(RefAssemblyIs(0, tseq("abcdefghijklm").size())));
}

TEST_F(discover_test, multiple) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd"), tseq("defgh"), tseq("deFgh"), tseq("ghij")});

  start();
  add_ref_asm(0, tseq("abcde"));
  add_ref_asm(tseq("abcde").size(), tseq("f"), [](const active_assembly* act) {
    EXPECT_THAT(act->a->rc_seqset_entries.ends(), ElementsAre(SeqsetEntryIs(tseq_rc("de"))));
    EXPECT_THAT(act->a->rc_seqset_entries.starts(), ElementsAre(SeqsetEntryIs(tseq_rc("def"))));
  });
  add_var_asm(
      tseq("abcde").size(), tseq("F"), tseq("abcdef").size(), [](const active_assembly* act) {
        EXPECT_THAT(act->a->rc_seqset_entries.ends(), ElementsAre(SeqsetEntryIs(tseq_rc("de"))));
        EXPECT_THAT(act->a->rc_seqset_entries.starts(), ElementsAre(SeqsetEntryIs(tseq_rc("deF"))));
      });
  add_ref_asm(tseq("abcdef").size(), tseq("ghij"), [](const active_assembly* act) {
    EXPECT_THAT(act->a->rc_seqset_entries.ends(),
                ElementsAre(SeqsetEntryIs(tseq_rc("def")), SeqsetEntryIs(tseq_rc("deF"))));
    EXPECT_THAT(act->a->rc_seqset_entries.starts(), ElementsAre(SeqsetEntryIs(tseq_rc("ghij"))));
  });
  flush();

  EXPECT_THAT(m_assemblies, SizeIs(4));
  EXPECT_THAT(m_traced, SizeIs(4));
}

TEST_F(discover_test, half_aligned) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  use_reads({tseq("abcd"), tseq("defgh"), tseq("de<gh"), tseq("de>gh"), tseq("ghij")});

  start();
  add_ref_asm(0, tseq("abcde"));
  add_ref_asm(tseq("abcde").size(), tseq("f"), [](const active_assembly* act) {
    EXPECT_THAT(act->a->rc_seqset_entries.ends(), ElementsAre(SeqsetEntryIs(tseq_rc("de"))));
    EXPECT_THAT(act->a->rc_seqset_entries.starts(), ElementsAre(SeqsetEntryIs(tseq_rc("def"))));
  });
  add_var_asm(
      tseq("abcde").size(), tseq(">"), optional_aoffset::none, [](const active_assembly* act) {
        EXPECT_THAT(act->a->rc_seqset_entries.ends(), ElementsAre(SeqsetEntryIs(tseq_rc("de"))));
        EXPECT_THAT(act->a->rc_seqset_entries.starts(), ElementsAre(SeqsetEntryIs(tseq_rc("de>"))));
      });
  add_var_asm(
      optional_aoffset::none, tseq("<"), tseq("abcdef").size(), [](const active_assembly* act) {
        EXPECT_THAT(act->a->rc_seqset_entries.ends(), ElementsAre(SeqsetEntryIs(tseq_rc(""))));
        EXPECT_THAT(act->a->rc_seqset_entries.starts(), ElementsAre(SeqsetEntryIs(tseq_rc("<"))));
      });
  add_ref_asm(tseq("abcdef").size(), tseq("ghij"), [](const active_assembly* act) {
    EXPECT_THAT(act->a->rc_seqset_entries.ends(),
                ElementsAre(SeqsetEntryIs(tseq_rc("def")), SeqsetEntryIs(tseq_rc("<"))));
    EXPECT_THAT(act->a->rc_seqset_entries.starts(), ElementsAre(SeqsetEntryIs(tseq_rc("ghij"))));
  });
  flush();

  EXPECT_THAT(m_assemblies, SizeIs(5));
  EXPECT_THAT(m_traced, SizeIs(5));
}

}  // namespace variants

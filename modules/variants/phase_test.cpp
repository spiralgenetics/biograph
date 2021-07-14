#include "modules/variants/phase.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class join_phases_test : public assemble_test {
 public:
  void add(const std::initializer_list<const char*>& phase_ids, aoffset_t left_offset,
           aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->tags.insert("phase_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = var_seq(phase_ids);
    a->phase_ids.insert(phase_ids.begin(), phase_ids.end());
    m_pipeline->add(std::move(a));
  }

  dna_sequence var_seq(const std::initializer_list<const char*>& phase_ids) {
    phase_set pids(phase_ids.begin(), phase_ids.end());
    return tseq(pids.to_string());
  }

  dna_sequence ref_seq(aoffset_t left_offset, aoffset_t right_offset) {
    return dna_sequence(
        m_scaffold.subscaffold(left_offset, right_offset - left_offset).get_simple());
  }

  join_phases_test() { m_chunk_bases = tseq("0").size(); }

  void start_join_phases(aoffset_t max_phase_len, aoffset_t max_phase_asm_len, size_t ref_pad) {
    join_phases::g_check_invariants = true;
    use_ref_parts({{0, tseq("0123456789")}});

    pipeline_step_t join =
        make_unique<join_phases>(max_phase_len, max_phase_asm_len, test_output());
    m_pipeline = make_unique<add_ref>(m_options, ref_pad, false /* whole ref */,
                                      0 /* length limit */, std::move(join));
  }

  void flush() {
    m_pipeline.reset();
    expect_sorted(assembly::left_offset_less_than);
  }

 protected:
  size_t m_chunk_bases = 0;
  pipeline_step_t m_pipeline;
};

std::vector<dna_sequence> get_sub_seqs(const assembly& a) {
  std::vector<dna_sequence> seqs;
  for (const auto& sub : a.sub_assemblies) {
    seqs.push_back((*sub)->seq);
  }
  return seqs;
}

/*Matcher<assembly> JoinedPhaseIs(std::vector<std::string> phase_ids, aoffset_t left_offset,
                                aoffset_t right_offset, std::vector<dna_sequence> sub_seqs) {
  return AllOf(Field(&assembly::phase_ids, ElementsAreArray(phase_ids)),
               Field(&assembly::left_offset, Eq(left_offset)),
               Field(&assembly::right_offset, Eq(right_offset)),
               ResultOf(get_sub_seqs, ElementsAreArray(sub_seqs)));
               }*/

MATCHER_P4(JoinedPhaseIsInternal, expected_phase_ids, expected_left_offset, expected_right_offset,
           expected_sub_seqs, "") {
  bool ok = true;

  if (expected_phase_ids != arg.phase_ids) {
    (*result_listener) << " where the phase ids " << arg.phase_ids << " is not the expected "
                       << expected_phase_ids;
    ok = false;
  }
  if (expected_left_offset != arg.left_offset) {
    (*result_listener) << " where the left offset " << arg.left_offset << " is not the expected "
                       << expected_left_offset;
    ok = false;
  }
  if (expected_right_offset != arg.right_offset) {
    (*result_listener) << " where the right offset " << arg.right_offset << " is not the expected "
                       << expected_right_offset;
    ok = false;
  }
  if (!ok) {
    *result_listener << "\n";
  }
  if (!ExplainMatchResult(ContainerEq(expected_sub_seqs), get_sub_seqs(arg), result_listener)) {
    ok = false;
  }
  return ok;
}

Matcher<assembly> JoinedPhaseIs(phase_set phase_ids, aoffset_t left_offset, aoffset_t right_offset,
                                std::vector<dna_sequence> sub_seqs) {
  return JoinedPhaseIsInternal(phase_ids, left_offset, right_offset, sub_seqs);
}

MATCHER_P(PhasedRefIs, seq_arg, "") {
  dna_sequence seq(seq_arg);
  bool ok = true;

  if (!arg.matches_reference) {
    ok = false;
    (*result_listener) << " where the assembly is not a reference assembly";
  }

  if (arg.seq != seq) {
    ok = false;
    (*result_listener) << " where the subassembly sequence " << arg.seq << " doesn't match " << seq;
  }

  if (arg.sub_assemblies.size() != 1) {
    ok = false;
    (*result_listener) << " which does not have exactly one subassembly";
  } else {
    const auto& sub = *arg.sub_assemblies[0];
    if (!sub->matches_reference) {
      ok = false;
      (*result_listener) << " where the subassembly does not match reference";
    }

    if (sub->seq != seq) {
      ok = false;
      (*result_listener) << " where the subassembly sequence " << sub->seq << " doesn't match "
                         << seq;
    }
  }

  return ok;
}

TEST_F(join_phases_test, trivial) {
  start_join_phases(1000, 1000, 0);
  add({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
  flush();

  EXPECT_THAT(
      m_non_ref_assemblies,
      ElementsAre(JoinedPhaseIs({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases, {var_seq({"a"})})));
  EXPECT_THAT(m_ref_assemblies, ElementsAre(PhasedRefIs(tseq("1"))));
}

TEST_F(join_phases_test, simple_join) {
  start_join_phases(m_chunk_bases, 1000, m_chunk_bases);
  add({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
  add({"a"}, 3 * m_chunk_bases, 4 * m_chunk_bases);
  flush();

  EXPECT_THAT(
      m_non_ref_assemblies,
      ElementsAre(JoinedPhaseIs(
          {"a"}, 1 * m_chunk_bases, 4 * m_chunk_bases,
          {var_seq({"a"}), ref_seq(2 * m_chunk_bases, 3 * m_chunk_bases), var_seq({"a"})})));
  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(PhasedRefIs(tseq("0")),  //
                          PhasedRefIs(tseq("1")),  //
                          PhasedRefIs(tseq("2")),  //
                          PhasedRefIs(tseq("3")),  //
                          PhasedRefIs(tseq("4"))));
}

TEST_F(join_phases_test, max_phase_dist) {
  start_join_phases(m_chunk_bases - 1, 1000, m_chunk_bases);
  add({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
  add({"a"}, 3 * m_chunk_bases, 4 * m_chunk_bases);
  flush();

  EXPECT_THAT(
      m_non_ref_assemblies,
      ElementsAre(JoinedPhaseIs({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases, {var_seq({"a"})}),
                  JoinedPhaseIs({"a"}, 3 * m_chunk_bases, 4 * m_chunk_bases, {var_seq({"a"})})));
  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(PhasedRefIs(tseq("0")),  //
                          PhasedRefIs(tseq("1")),  //
                          PhasedRefIs(tseq("2")),  //
                          PhasedRefIs(tseq("3")),  //
                          PhasedRefIs(tseq("4"))));
}

TEST_F(join_phases_test, adjacent) {
  start_join_phases(1000, 1000, m_chunk_bases);
  add({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
  add({"a"}, 2 * m_chunk_bases, 3 * m_chunk_bases);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(JoinedPhaseIs({"a"}, 1 * m_chunk_bases, 3 * m_chunk_bases,
                                        {var_seq({"a"}), var_seq({"a"})})));
  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(PhasedRefIs(tseq("0")),  //
                          PhasedRefIs(tseq("1")),  //
                          PhasedRefIs(tseq("2")),  //
                          PhasedRefIs(tseq("3"))));
}

TEST_F(join_phases_test, conflict) {
  start_join_phases(1000, 1000, m_chunk_bases);
  EXPECT_THROW(
      {
        add({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
        add({"a"}, 2 * m_chunk_bases - 1, 3 * m_chunk_bases);
        // Have to cause it to flush here, since exceptions in destructors
        // are not propagated:
        add({"b"}, 10 * m_chunk_bases, 11 * m_chunk_bases);
        flush();
      },
      io_exception);
}

TEST_F(join_phases_test, complex_conflict) {
  start_join_phases(1000, 1000, m_chunk_bases);
  EXPECT_THROW(
      {
        add({"a", "b"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
        add({"a"}, 2 * m_chunk_bases - 1, 3 * m_chunk_bases);
        add({"b"}, 3 * m_chunk_bases, 4 * m_chunk_bases);
        flush();
      },
      io_exception);
}

TEST_F(join_phases_test, split) {
  start_join_phases(1000, 1000, m_chunk_bases);
  add({"a", "b"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
  add({"b"}, 3 * m_chunk_bases, 4 * m_chunk_bases);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(  //

                  JoinedPhaseIs({"b"}, 1 * m_chunk_bases, 4 * m_chunk_bases,
                                {                                                //
                                 var_seq({"a", "b"}),                            //
                                 ref_seq(2 * m_chunk_bases, 3 * m_chunk_bases),  //
                                 var_seq({"b"})}),
                  JoinedPhaseIs({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases,
                                {//
                                 var_seq({"a", "b"})})));
  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(PhasedRefIs(tseq("0")),  //
                          PhasedRefIs(tseq("1")),  //
                          PhasedRefIs(tseq("2")),  //
                          PhasedRefIs(tseq("3")),  //
                          PhasedRefIs(tseq("4"))));
}

TEST_F(join_phases_test, complex_split) {
  start_join_phases(1000, 1000, m_chunk_bases);
  add({"a", "b", "c"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
  add({"b", "c"}, 3 * m_chunk_bases, 4 * m_chunk_bases);
  add({"a", "d"}, 5 * m_chunk_bases, 6 * m_chunk_bases);
  flush();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(  //
                  JoinedPhaseIs({"a"}, 1 * m_chunk_bases, 6 * m_chunk_bases,
                                {                                                //
                                 var_seq({"a", "b", "c"}),                       //
                                 ref_seq(2 * m_chunk_bases, 3 * m_chunk_bases),  //
                                 ref_seq(3 * m_chunk_bases, 4 * m_chunk_bases),  //
                                 ref_seq(4 * m_chunk_bases, 5 * m_chunk_bases),  //
                                 var_seq({"a", "d"})}),
                  JoinedPhaseIs({"b", "c"}, 1 * m_chunk_bases, 4 * m_chunk_bases,
                                {                                                //
                                 var_seq({"a", "b", "c"}),                       //
                                 ref_seq(2 * m_chunk_bases, 3 * m_chunk_bases),  //
                                 var_seq({"b", "c"})}),
                  JoinedPhaseIs({"d"}, 5 * m_chunk_bases, 6 * m_chunk_bases,
                                {//
                                 var_seq({"a", "d"})})));
  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(PhasedRefIs(tseq("0")),  //
                          PhasedRefIs(tseq("1")),  //
                          PhasedRefIs(tseq("2")),  //
                          PhasedRefIs(tseq("3")),  //
                          PhasedRefIs(tseq("4")),  //
                          PhasedRefIs(tseq("5")),  //
                          PhasedRefIs(tseq("6"))));
}

class split_phases_test : public assemble_test {
 public:
  void add(const std::initializer_list<const char*>& phase_ids, aoffset_t left_offset,
           aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    m_expected_asm_ids.insert(a->assembly_id);
    a->tags.insert("phase_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = tseq(PrintToString(phase_ids));
    a->phase_ids.insert(phase_ids.begin(), phase_ids.end());
    m_pipeline->add(std::move(a));
  }

  void start_split_phases(aoffset_t max_phase_len, size_t ref_pad) {
    join_phases::g_check_invariants = true;
    use_ref_parts({{0, tseq("0123456789")}});

    pipeline_step_t split = make_unique<split_phases>(test_output());
    pipeline_step_t join = make_unique<join_phases>(max_phase_len, 1000, std::move(split));
    m_pipeline = make_unique<add_ref>(m_options, ref_pad, false /* whole ref */,
                                      0 /* length limit */, std::move(join));
  }

  void flush() {
    m_pipeline.reset();
    expect_sorted(assembly::left_offset_less_than);

    for (const auto& a : m_non_ref_assemblies) {
      m_asm_ids.insert(a.assembly_id);
    }

    EXPECT_THAT(m_asm_ids, ContainerEq(m_expected_asm_ids));
  }

  split_phases_test() { m_chunk_bases = tseq("0").size(); }

 protected:
  pipeline_step_t m_pipeline;
  size_t m_chunk_bases = 0;

  std::multiset<size_t> m_expected_asm_ids;
  std::multiset<size_t> m_asm_ids;
};

TEST_F(split_phases_test, split) {
  start_split_phases(1000, m_chunk_bases);
  add({"a", "b", "c"}, 1 * m_chunk_bases, 2 * m_chunk_bases);
  add({"b", "c"}, 3 * m_chunk_bases, 4 * m_chunk_bases);
  add({"a", "d"}, 5 * m_chunk_bases, 6 * m_chunk_bases);
  flush();
}

class resolve_conflicts_test : public assemble_test {
 public:
  resolve_conflicts_test() {
    m_do_resolve = [this](const assembly_ptr& a, const assembly_ptr& b, const phase_set& ids) {
      m_conflicts[a->assembly_id]++;
      m_conflicts[b->assembly_id]++;
      for (const auto& id : ids) {
        if (a->assembly_id < b->assembly_id) {
          a->phase_ids.erase(id);
        } else {
          b->phase_ids.erase(id);
        }
      }
    };
    m_chunk_bases = tseq("0").size();
  }
  void add(const std::initializer_list<const char*>& phase_ids, aoffset_t left_offset,
           aoffset_t right_offset, size_t asm_id) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = asm_id;
    m_expected_asm_ids.insert(a->assembly_id);
    a->tags.insert("phase_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = tseq(PrintToString(phase_ids));
    a->phase_ids.insert(phase_ids.begin(), phase_ids.end());
    m_pipeline->add(std::move(a));
  }

  void start_resolve_conflicts() {
    m_pipeline = ::make_unique<resolve_phase_conflicts>(m_do_resolve, test_output());
  }

  void flush() {
    m_pipeline.reset();
    expect_sorted(assembly::left_offset_less_than);

    for (const auto& a : m_non_ref_assemblies) {
      m_asm_ids.insert(a.assembly_id);
    }

    EXPECT_THAT(m_asm_ids, ContainerEq(m_expected_asm_ids));
  }

 protected:
  pipeline_step_t m_pipeline;

  resolve_phase_conflicts::resolve_conflict_func_t m_do_resolve;
  std::map<size_t /* assembly id */, size_t /* conflict count */> m_conflicts;

  std::multiset<size_t> m_expected_asm_ids;
  std::multiset<size_t> m_asm_ids;
  size_t m_chunk_bases = 0;
};

TEST_F(resolve_conflicts_test, simple_conflict) {
  start_resolve_conflicts();
  add({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases, 1);
  add({"a"}, 2 * m_chunk_bases - 1, 3 * m_chunk_bases, 2);
  flush();

  EXPECT_THAT(m_conflicts, ElementsAre(Pair(1, 1), Pair(2, 1)));
}

TEST_F(resolve_conflicts_test, complex_conflict) {
  start_resolve_conflicts();
  add({"a", "b"}, 1 * m_chunk_bases, 2 * m_chunk_bases, 1);
  add({"a"}, 2 * m_chunk_bases - 1, 3 * m_chunk_bases, 2);
  add({"b"}, 3 * m_chunk_bases, 4 * m_chunk_bases, 3);
  flush();

  EXPECT_THAT(m_conflicts, ElementsAre(Pair(1, 1), Pair(2, 1)));
}

TEST_F(resolve_conflicts_test, fails_to_resolve) {
  m_do_resolve = [](const assembly_ptr& a, const assembly_ptr& b, const phase_set& phase_ids) {
    // Don't do anything to resolve the conflict.
  };
  start_resolve_conflicts();
  EXPECT_THROW(
      {
        add({"a"}, 1 * m_chunk_bases, 2 * m_chunk_bases, 1);
        add({"a"}, 2 * m_chunk_bases - 1, 3 * m_chunk_bases, 2);
        add({"b"}, 10 * m_chunk_bases - 1, 11 * m_chunk_bases, 3);
        flush();
      },
      std::runtime_error);
}

}  // namespace variants

#include "modules/variants/pipeline.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/trace_ref.h"
#include "modules/variants/align.h"
#include "modules/variants/ploid_limit.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class pipeline_test
    : public assemble_test,
      public WithParamInterface<std::pair<bool /* rev_comp */, bool /* bidir tracer */>> {
 public:
  pipeline_test() { std::tie(m_rev_comp, m_bidir_tracer)  = GetParam();  }
  static std::string get_ref_scaffold_name(const dna_sequence& seq,
                                           size_t expected_position = 0) {
    const auto& flat_ref = m_ref->get_flat_ref();
    const flat_ref::index_t& index = flat_ref.get_index();
    for (const auto& extent : index.extents) {
      dna_slice ext_slice(m_ref->get_dna(extent.flat), extent.size);

      std::cout << "Found ext slice: " << ext_slice << "\n";
      std::cout.flush();
      if (ext_slice == seq) {
        seq_position pos = m_ref->get_seq_position(extent.flat);
        CHECK_EQ(expected_position, pos.position);
        return m_ref->get_assembly().scaffold_order[pos.scaffold_id];
      }
    }

    LOG(FATAL) << "No scaffold found for sequence " << seq;
    return "";
  }

  static void make_reference() {
    m_ref_contents =                                       //
        tseq("abc").as_string() + std::string(100, 'N') +  //
        tseq("defg").as_string() + "G" +                   //
        tseq("hijkl").as_string() + "GG" +                 //
        tseq("mnopqrstuvwxyz").as_string() +               //
        // Repetitive region
        tseq("ABCD").as_string() + "A" + tseq("EFGH").as_string() +  //
        "T" + tseq("EFGHIJKL").as_string();

    m_ref = create_reference_str({m_ref_contents});
    m_def_start = tseq("abc").size() + 100;
    m_hij_start = m_def_start + (tseq("defg") + dna_G).size();
    m_mno_start = m_hij_start + (tseq("hijkl") + dna_G + dna_G).size();
    m_ABC_start = m_mno_start + tseq("mnopqrstuvwxyz").size();
    m_scaffold_name = get_ref_scaffold_name(tseq("abc"));
  }

  void run_pipeline() {
    if (!m_ref) {
      make_reference();
    }
    m_scaffold = trace_ref::ref_to_scaffold(m_ref.get(), m_scaffold_name);
    CHECK_EQ(m_scaffold.as_string(), m_ref_contents);
    CHECK_EQ(m_scaffold.subscaffold_str(m_def_start, tseq("def").size()), tseq("def").as_string());
    CHECK_EQ(m_scaffold.subscaffold_str(m_hij_start, tseq("hij").size()), tseq("hij").as_string());
    CHECK_EQ(m_scaffold.subscaffold_str(m_mno_start, tseq("mno").size()), tseq("mno").as_string());
    CHECK_EQ(m_scaffold.subscaffold_str(m_ABC_start, tseq("ABC").size()), tseq("ABC").as_string());

    m_rmap.emplace(m_seqset.get(), m_ref.get());
    m_rmap->build();

    m_options.ref = m_ref.get();
    m_options.rmap = &m_rmap.get();
    m_options.use_bidir_tracer = m_bidir_tracer;
    m_options.debug_paths = [&](const std::string& dot_contents) {
      std::string filename = "/tmp/path-debug.dot";
      static std::atomic<size_t> next_debug{0};
      filename += ".";
      filename += std::to_string(next_debug.fetch_add(1));
      std::cout << "Writing path debug to " << filename << "\n";
      std::ofstream path_debug(filename);
      path_debug << dot_contents;
    };
    // TODO(nils): Remove this option and rework the test not to need it.
    m_options.trace_reference_assemblies = true;

    m_options.min_overlap = 2 * k_dna_test_sequence_length;
    m_pipeline.emplace(m_options, test_output());
    m_pipeline->add_step<aligner>(m_options);
    m_pipeline->add_step<ploid_limiter>(m_options);
    // TODO(nils): Rework test to not require the align splitter.
    m_pipeline->add_step<align_splitter>();

    m_options.scaffold = nullptr;
    test_scaffold_pipeline sp(m_scaffold_name, &m_pipeline.get());
    m_trace.emplace(m_options, &sp);
    if (!m_bidir_tracer) {
      if (m_rev_comp) {
        m_options.skip_push_trace_fwd = true;
      } else {
        m_options.skip_push_trace_rev = true;
      }
    }
    m_trace->add_scaffold(m_scaffold_name);
    auto st = m_trace->assemble();
    std::cout << "Assembly complete; stats: " << st << "\n";
    m_trace.reset();
    m_pipeline.reset();
  }

 protected:
  boost::optional<assemble_pipeline> m_pipeline;
  boost::optional<ref_map> m_rmap;
  boost::optional<trace_ref> m_trace;
  bool m_rev_comp;
  bool m_bidir_tracer;

  static std::unique_ptr<reference> m_ref;
  static std::string m_ref_contents;
  static aoffset_t m_def_start;
  static aoffset_t m_hij_start;
  static aoffset_t m_mno_start;
  static aoffset_t m_ABC_start;
  static std::string m_scaffold_name;
};

std::unique_ptr<reference> pipeline_test::m_ref;
std::string pipeline_test::m_ref_contents;
aoffset_t pipeline_test::m_def_start;
aoffset_t pipeline_test::m_hij_start;
aoffset_t pipeline_test::m_mno_start;
aoffset_t pipeline_test::m_ABC_start;
std::string pipeline_test::m_scaffold_name;

TEST_P(pipeline_test, ref_only) {
  use_reads({tseq("mnopq"), tseq("opqrst"), tseq("rstuv"), tseq("uvwxy")});
  run_pipeline();

  EXPECT_THAT(m_non_ref_assemblies, IsEmpty());
}

TEST_P(pipeline_test, homozygous_snp) {
  g_trace_all_assemblies = true;
  use_reads({tseq("defg"), tseq("fg") + dna_C + tseq("hi"), tseq("hijk")});
  run_pipeline();

  EXPECT_THAT(
      m_non_ref_assemblies,
      ElementsAre(AssemblyIs(m_def_start + tseq("defg").size(), dna_C,
                             m_def_start + (tseq("defg") + dna_C).size())));
}

TEST_P(pipeline_test, hetrozygous_snp) {
  use_reads({tseq("defg"),
             // Variant:
             tseq("fg") + dna_A + tseq("hi"),
             // Reference:
             tseq("fg") + dna_G + tseq("hi"), tseq("hijk")});
  run_pipeline();

  EXPECT_THAT(
      m_non_ref_assemblies,
      ElementsAre(AssemblyIs(m_def_start + tseq("defg").size(), dna_A,
                             m_def_start + (tseq("defg") + dna_A).size())));
}

TEST_P(pipeline_test, compound_hetrozygous_snp) {
  use_reads({tseq("defg"),
             // Variant 1:
             tseq("fg") + dna_A + tseq("hi"),
             // Variant 2:
             tseq("fg") + dna_T + tseq("hi"), tseq("hijk")});
    run_pipeline();


  EXPECT_THAT(m_non_ref_assemblies,
              UnorderedElementsAre(
                  AssemblyIs(m_def_start + tseq("defg").size(), dna_A,
                             m_def_start + (tseq("defg") + dna_G).size()),
                  AssemblyIs(m_def_start + tseq("defg").size(), dna_T,
                             m_def_start + (tseq("defg") + dna_G).size())));
}

TEST_P(pipeline_test, interscaffold_delete) {
  use_reads({tseq("abc"), tseq("bcd"), tseq("cde"), tseq("def"), tseq("efg")});
  run_pipeline();

  EXPECT_THAT(
      m_non_ref_assemblies,
      ElementsAre(AssemblyIs(tseq("abc").size(), dna_sequence(), m_def_start)));
}

TEST_P(pipeline_test, interscaffold_insert) {
  use_reads({tseq("abc"), tseq("bc") + dna_C + tseq("d"),
             tseq("c") + dna_C + tseq("de"), tseq("def"), tseq("efg")});
  run_pipeline();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AssemblyIs(tseq("abc").size(), dna_C, m_def_start)));
}

// Test the case where we could call the same thing more than one
// different way; we want to choose the simpler way.
//
// This was an attempt to reproduce DEV-407, but upon further
// investigation it seems like a different problem.
//
// TODO(nils): Make this test pass.
TEST_P(pipeline_test, DISABLED_self_similar_snp) {
  // Reference is tseq("ABCD") + dna_A + tseq("EFGH") + dna_T +
  // tseq("EFGH") Should be called as a SNP on the first dna_A instead
  // of a deletion of dna_A + tseq("EFGH").
  use_reads({tseq("ABCD"), tseq("BCD") + dna_T + tseq("E"),
          tseq("CD") + dna_T + tseq("EF"), tseq("D") + dna_T + tseq("EFG"),
          dna_T + tseq("EFGH")});
  run_pipeline();

  EXPECT_THAT(m_non_ref_assemblies,
              ElementsAre(AssemblyIs(m_ABC_start + tseq("ABCD").size(),
                                     dna_T,
                                     m_ABC_start + tseq("ABCD").size() + 1)));
}

INSTANTIATE_TEST_CASE_P(fwd_pipeline_test, pipeline_test,
                        ::testing::Values(std::make_pair(false /* not rev_comp */,
                                                         false /* not bidir */)));
INSTANTIATE_TEST_CASE_P(rev_pipeline_test, pipeline_test,
                        ::testing::Values(std::make_pair(true /*  rev_comp */,
                                                         false /* not bidir */)));
INSTANTIATE_TEST_CASE_P(bidir_pipeline_test, pipeline_test,
                        ::testing::Values(std::make_pair(false, true /* bidir */)));

}  // namespace variants

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/ref_map.h"
#include "modules/variants/trace_ref.h"

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

class trace_ref_test : public Test, public pipeline_interface {
 protected:
  static void make_reference() {
    auto num = tseq("0123456789");
    auto alpha1 = tseq("abcdefg");
    auto alpha2 = tseq("hijklm");
    auto alpha3 = tseq("nopqrstuvw");
    // Hvae to have a fair number of 'N's in order to actually get
    // multiple extents; see fast_ref_importer::add_base.
    auto alpha = alpha1.as_string() + std::string(alpha2.size(), 'N') +
                 alpha3.as_string();
    m_ref = create_reference_str({alpha, num.as_string()});
    m_num_scaffold_name = get_ref_scaffold_name(num);
    m_alpha_scaffold_name = get_ref_scaffold_name(alpha1);
    CHECK_EQ(m_alpha_scaffold_name,
             get_ref_scaffold_name(alpha3, alpha1.size() + alpha2.size()));

    CHECK_EQ(alpha1, get_ref_part_seq(m_alpha_scaffold_name, 0, alpha1.size()));
    CHECK_EQ(alpha3,
             get_ref_part_seq(m_alpha_scaffold_name,
                              alpha1.size() + alpha2.size(), alpha3.size()));
    CHECK_EQ(num, get_ref_part_seq(m_num_scaffold_name, 0, num.size()));
  }

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

  trace_ref_test() {
    if (!m_ref) {
      make_reference();
    }
    CHECK(m_ref);

    m_opts.min_overlap = k_dna_test_sequence_length * 2;
    // TODO(nils): Remove these option and rework the test not to need them.
    m_opts.use_bidir_tracer = false;
    m_opts.trace_reference_assemblies = true;
    m_opts.pop_trace_anchor_drop = false;
    m_opts.report_half_aligned_func = [&](
        const half_aligned_assembly& ha) {
      std::cout << "Got half-aligned: " << ha << "\n";
    };
    g_trace_all_assemblies = true;
  }

  pipeline_step_t make_parallel_input() override {
    return make_unique<assemble_lambda_output>(
        [&](assembly_ptr a) {
          static std::mutex mu;
          std::lock_guard<std::mutex> l(mu);
          std::cout << "Got assembly: " << *a << "\n";
          m_assemblies.push_back(*a);

          ASSERT_GT(a->seq.size(), a->left_anchor_len);
          ASSERT_GT(a->seq.size(), a->right_anchor_len);
          EXPECT_EQ(a->seq.subseq(0, a->left_anchor_len),
                    get_ref_part_seq(m_cur_scaffold_name, a->left_offset,
                                     a->left_anchor_len))
              << *a;
          EXPECT_EQ(
              a->seq.subseq(aoffset_t(a->seq.size()) - a->right_anchor_len,
                            a->right_anchor_len),
              get_ref_part_seq(m_cur_scaffold_name,
                               a->right_offset - a->right_anchor_len,
                               a->right_anchor_len))
              << *a;
        },
        "raw_assemblies");
  }

  void do_assemble(std::string scaffold_name) {
    m_cur_scaffold_name = scaffold_name;
    m_seqset = seqset_for_reads(m_all_reads);
    m_readmap = readmap_for_reads(m_seqset, m_paired_reads, m_reads);

    if (!m_ref) {
      make_reference();
    }
    CHECK(m_ref);
    m_rmap.emplace(m_seqset.get(), m_ref.get());
    m_rmap->build();

    m_assemblies.clear();
    m_opts.seqset = m_seqset.get();
    m_opts.readmap = m_readmap.get();
    m_opts.ref = m_ref.get();
    m_opts.rmap = &m_rmap.get();

    test_scaffold_pipeline p(scaffold_name, this);
    m_trace.emplace(m_opts, &p);
    m_trace->add_scaffold(scaffold_name);
    auto st = m_trace->assemble();
    std::cout << "Assemble stats: " << st << "\n";
    m_trace.reset();
  }

  static dna_sequence get_ref_part_seq(std::string scaffold_name,
                                       aoffset_t offset, int len) {
    std::cout << "Getting part seq for " << offset << ", len " << len << "\n";
    const auto& refasm = m_ref->get_assembly();
    size_t flat_start = refasm.flatten(scaffold_name, offset);
    size_t flat_end = refasm.flatten(scaffold_name, offset + len);
    return dna_sequence(m_ref->get_dna(flat_start), m_ref->get_dna(flat_end));
  }

  void add_reads(const std::vector<dna_sequence>& reads) {
    for (const auto& read : reads) {
      m_reads.push_back(read);
      m_all_reads.push_back(read);
    }
  }
  void add_paired_reads(const std::vector<std::pair<dna_sequence, dna_sequence>>& reads) {
    add_paired_reads(reads, {});
  }
  void add_paired_reads(const std::vector<std::pair<dna_sequence, dna_sequence>>& reads, const std::vector<dna_sequence>& unpaired) {
    for (const auto& read : reads) {
      m_paired_reads.push_back(read);
      m_all_reads.push_back(read.first);
      m_all_reads.push_back(read.second);
    }
    for (const auto& read : unpaired) {
      m_reads.push_back(read);
      m_all_reads.push_back(read);
    }
  }

  assemble_options m_opts;
  std::string m_cur_scaffold_name;
  static std::string m_alpha_scaffold_name, m_num_scaffold_name;
  pipeline_step_t m_output_f;
  boost::optional<ref_map> m_rmap;
  static std::unique_ptr<reference> m_ref;
  std::shared_ptr<seqset> m_seqset;
  std::unique_ptr<readmap> m_readmap;
  std::vector<assembly> m_assemblies;
  boost::optional<trace_ref> m_trace;
  std::vector<dna_sequence> m_reads;
  std::vector<dna_sequence> m_all_reads;
  std::vector<std::pair<dna_sequence, dna_sequence>> m_paired_reads;
  scaffold m_scaffold;
};

class trace_ref_push_test : public trace_ref_test, public WithParamInterface<bool /* rev_comp */> {
 public:
  trace_ref_push_test() {
    m_rev_comp = GetParam();
    if (m_rev_comp) {
      m_opts.skip_push_trace_fwd = true;
    } else {
      m_opts.skip_push_trace_rev = true;
    }
  }

 protected:
  bool m_rev_comp;
};

std::unique_ptr<reference> trace_ref_test::m_ref;
std::string trace_ref_test::m_alpha_scaffold_name,
    trace_ref_test::m_num_scaffold_name;

TEST_P(trace_ref_push_test, all_reference) {
  add_reads({tseq("012345"), tseq("34567")});
  do_assemble(m_num_scaffold_name);
  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(RefTrAssemblyIs(tseq("01234567"))));
}

TEST_P(trace_ref_push_test, spans_extent) {
  add_reads({tseq("bcdef"), tseq("efghijklmnop"), tseq("nopqrstuvw")});
  do_assemble(m_alpha_scaffold_name);
  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(TrAssemblyIs(tseq("bcdef"), tseq("ghijklm"),
                                                tseq("nopqrstuvw"))));
}


INSTANTIATE_TEST_CASE_P(fwd_trace_ref_push_test, trace_ref_push_test,
                        ::testing::Values(false /* not rev_comp */));
INSTANTIATE_TEST_CASE_P(rev_trace_ref_push_test, trace_ref_push_test,
                        ::testing::Values(true /*  rev_comp */));

class trace_ref_pop_test : public trace_ref_test, public WithParamInterface<bool /* rev_comp */> {
 public:
  trace_ref_pop_test() {
    CHECK(m_ref);

    m_rev_comp = GetParam();
    if (m_rev_comp) {
      m_opts.skip_pop_trace_fwd = true;
    } else {
      m_opts.skip_pop_trace_rev = true;
    }

    m_opts.pop_trace_anchor_drop = true;
    m_opts.forward_pairs_face_inward = false;
    m_opts.min_overlap = k_dna_test_sequence_length * 3;
    m_opts.min_pair_distance = 10;
    m_opts.max_pair_distance = 1000;
    m_opts.min_pop_overlap = k_dna_test_sequence_length * 2;
  }

 protected:
  bool m_rev_comp;
};

TEST_P(trace_ref_pop_test, simple) {
  add_paired_reads({//
                    // full reads before and after variant, which the push tracer needs to anchor.
                    {tseq_rc("nopq"), tseq("tuvw")},
                    // First side:
                    {tseq_rc("abcde"), tseq("opq") + dna_G + tseq("R")},
                    {tseq_rc("abcde"), tseq("pq") + dna_G + tseq("RS")},

                    // Second side:
                    {tseq_rc("abcde"), tseq("RS") + dna_G + tseq("tu")},
                    {tseq_rc("abcde"), tseq("S") + dna_G + tseq("tuv")}});

  do_assemble(m_alpha_scaffold_name);
  EXPECT_THAT(m_assemblies,
              ElementsAre(AnyOf(
                  AssemblyIs(tseq("abcdefghijklmnopq").size(),
                             dna_G + tseq("RS") + dna_G + tseq("tuvw"),
                             tseq("abcdefghijklmnopqrstuvw").size()),
                  AssemblyIs(tseq("abcdefghijklm").size(),
                             tseq("nopq") + dna_G + tseq("RS") + dna_G,
                             tseq("abcdefghijklmnopqrs").size()))));
}

TEST_P(trace_ref_pop_test, additional_inside_read) {
  add_paired_reads({// Has a mate which supplies tseq("0123") to the pop tracer as a potential read.
                    // Anchor for push tracer:
                    {tseq_rc("nopq"), tseq("tuvw")},
                    // Supply pair information so it's aware of tseq("0123"):
                    {tseq_rc("nopq"), tseq("0123")},
                    // First side, with inside:
                    {tseq_rc("abcd"), tseq("opq") + dna_G + tseq("0")},
                    {tseq_rc("abcd"), tseq("pq") + dna_G + tseq("01")},

                    // Second side:
                    {tseq_rc("abcd"), tseq("23") + dna_G + tseq("tu")},
                    {tseq_rc("abcd"), tseq("3") + dna_G + tseq("tuv")}});

  do_assemble(m_alpha_scaffold_name);
  EXPECT_THAT(
      m_assemblies,
      UnorderedElementsAre(AnyOf(
          AssemblyIs(tseq("abcdefghijklmnopq").size(), dna_G + tseq("0123") + dna_G + tseq("tuvw"),
                     tseq("abcdefghijklmnopqrstuvw").size()),
          AssemblyIs(tseq("abcdefghijklm").size(), tseq("nopq") + dna_G + tseq("0123") + dna_G,
                     tseq("abcdefghijklmnopqrs").size()))));
}

INSTANTIATE_TEST_CASE_P(fwd_trace_ref_pop_test, trace_ref_pop_test,
                        ::testing::Values(false /* not rev_comp */));
INSTANTIATE_TEST_CASE_P(rev_trace_ref_pop_test, trace_ref_pop_test,
                        ::testing::Values(true /*  rev_comp */));


}  // namespace variants

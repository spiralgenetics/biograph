#pragma once

#include "modules/bio_base/biograph_dir.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/variants/assemble.h"
#include "modules/variants/assembly_dot.h"
#include "modules/variants/pipeline.h"
#include "modules/variants/scaffold.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace variants {

// Utility for running assemble tests on a large sample with real world data.
//
// Biographs referenced may be searched for in a number of different
// directories, so that faster local copies can be used.
//
// /scratch is assumed to be writable as the refmap will be stored there.
class big_assemble_test : public ::testing::Test, public pipeline_interface {
 protected:
  big_assemble_test() {
    disable_test_sequence_expansion();
    assemble_pipeline_interface::global_set_verify_order(true);
    std::cout.setf(std::ios::unitbuf);
    std::cerr.setf(std::ios::unitbuf);
  }
  void use_biograph(std::string bg_dir);

  void call_at(const std::string& scaffold_name, std::string one_based_pos,
               aoffset_t read_around_before, aoffset_t read_around_after);

  void call_region(const std::string& scaffold_name, aoffset_t start, aoffset_t limit);
  pipeline_step_t test_output();

  pipeline_step_t make_parallel_input() override;

  // Set tracing options such that we try extra hard.
  void set_thorough_trace_options();

  // Call at the given location and expect the given VCF entry.
  void run_vcf_test(const std::string& scaffold_name,
                    const std::string& vcf_start_offset /* 1-based offset */,
                    dna_sequence_matcher ref_bases, dna_sequence_matcher alt_1,
                    const std::string& gt_1);
  void run_vcf_test(const std::string& scaffold_name,
                    const std::string& vcf_start_offset /* 1-based offset */,
                    dna_sequence_matcher ref_bases, dna_sequence_matcher alt_1,
                    const std::string& gt_1, dna_sequence_matcher alt_2, const std::string& gt_2);

  // Print the given assemblies, but filter for ones where either the
  // reference or sequence is at least given the size.
  std::string dump_sv_assemblies(int min_size);

  dna_sequence get_ref_part_seq(aoffset_t offset, int len) {
    scaffold sub = m_scaffold.subscaffold(offset, len);
    CHECK(sub.is_simple());
    dna_slice slice = sub.get_simple();
    return dna_sequence(slice.begin(), slice.end());
  }

  void add_trace(aoffset_t left_offset, aoffset_t right_offset, dna_sequence seq);
  void add_trace(aoffset_t left_offset, aoffset_t right_offset, std::string seq) {
    add_trace(left_offset, right_offset, dna_sequence(seq));
  }

 private:
  void run_vcf_test_internal(const std::string& scaffold_name,
                             const std::string& vcf_start_offset /* 1-based offset */,
                             dna_sequence_matcher ref_bases, dna_sequence_matcher alt_1,
                             const std::string& gt_1, boost::optional<dna_sequence_matcher> alt_2,
                             boost::optional<std::string> gt_2);
  void select_scaffold(const std::string& scaffold_name);
  void call_region_internal(const std::string& scaffold_name, aoffset_t start, aoffset_t limit);
  static void init_search_path();
  static void open_biograph(std::string bg_path);

 protected:
  static std::vector<std::string> g_search_path;
  static boost::optional<reference> m_ref;
  static std::shared_ptr<seqset> m_seqset;
  static boost::optional<readmap> m_readmap;
  static std::string m_refmap_path;
  static boost::optional<ref_map> m_rmap;
  static std::string m_cur_biograph_dir;

  variants::scaffold m_scaffold;

  std::mutex m_mu;
  assemble_options m_options;
  size_t m_flat_call_pos = std::numeric_limits<size_t>::max();
  aoffset_t m_call_around_len = 5;
  aoffset_t m_call_pos;
  dna_const_iterator m_call_ref_it;

  aoffset_t m_interesting_left_offset = 0;
  aoffset_t m_interesting_right_offset = 0;

  boost::optional<assembly_dot> m_assembly_dot;
  boost::optional<assembly_dot> m_aligned_dot;

  std::unique_ptr<assemble_pipeline> m_pipeline;
  std::vector<assembly> m_assemblies;
  std::vector<assembly> m_var_assemblies;
  std::vector<assembly> m_ref_assemblies;
  std::vector<assembly> m_aligned;
  assemble_stats m_stats;
  bool m_trace_enabled = false;
  std::set<std::tuple<aoffset_t, aoffset_t, dna_sequence>> m_traces;
};

#define SCOPED_BIG_ASM_TEST() SCOPED_TRACE("Big assemble test")

MATCHER_P2(RefAt, expected_offset, expected_seq, "") {
  dna_sequence eseq(expected_seq);
  if (arg.aligned_variants.empty()) {
    // This version of this checker requires a reference-only assembly.
    if (!arg.matches_reference) {
      return false;
    }
    if (arg.right_offset < (expected_offset + aoffset_t(eseq.size())) ||
        arg.left_offset > expected_offset) {
      return false;
    }
    int left_bound = expected_offset - arg.left_offset;
    int right_bound = left_bound + eseq.size();
    dna_sequence assembly_seq = arg.seq.subseq(left_bound, right_bound - left_bound);
    if (assembly_seq == eseq) {
      return true;
    }
    *result_listener << "\nExpected: " << eseq << "\nIn assembly: " << arg.seq.subseq(0, left_bound)
                     << " [" << assembly_seq << "] "
                     << arg.seq.subseq(right_bound, arg.seq.size() - right_bound) << "\n";
    return false;
  } else {
    // Check the space between aligned variants to see if there's a
    // reference section that matches the expected sequence.
    aoffset_t offset_pos = arg.left_offset;
    aoffset_t seq_pos = 0;

    auto advance_and_check = [&](aoffset_t new_offset) -> bool {
      CHECK_GE(new_offset, offset_pos) << dump_assembly_and_vars(arg);
      aoffset_t adv = new_offset - offset_pos;
      CHECK_GE(adv, 0);
      aoffset_t seq_start = expected_offset - offset_pos + seq_pos;
      aoffset_t seq_end = seq_start + eseq.size();
      if (seq_start >= seq_pos && seq_end <= aoffset_t(arg.seq.size())) {
        dna_sequence sub = arg.seq.subseq(seq_start, eseq.size());
        if (sub == eseq) {
          return true;
        }

        *result_listener << "\nExpected: " << eseq << "\nIn: " << arg.seq.subseq(0, seq_start)
                         << " [ " << arg.seq.subseq(seq_start, eseq.size()) << " ] "
                         << arg.seq.subseq(seq_start + eseq.size(),
                                           arg.seq.size() - (seq_start + eseq.size()))
                         << " in assembly: " << dump_assembly_and_vars(arg) << "\n";
      }

      offset_pos += adv;
      seq_pos += adv;

      return false;
    };
    auto it = arg.aligned_variants.begin();
    while (it != arg.aligned_variants.end()) {
      if (advance_and_check(it->left_offset)) {
        return true;
      }
      seq_pos += it->seq.size();
      offset_pos += it->right_offset - it->left_offset;
      CHECK_EQ(offset_pos, it->right_offset);
      ++it;
    }
    if (advance_and_check(arg.right_offset)) {
      return true;
    }
    CHECK_EQ(offset_pos, arg.right_offset);
    return false;
  }
}

MATCHER_P3(VariantAt, offset, ref_size, expected_seq, "") {
  auto m = dna_sequence_matcher(expected_seq).matcher();
  if (m.Matches(arg.seq) && arg.left_offset == offset && arg.right_offset == offset + ref_size &&
      !arg.matches_reference) {
    return true;
  }
  for (const auto& v : arg.aligned_variants) {
    if (v.left_offset == offset && v.right_offset == offset + ref_size && m.Matches(v.seq)) {
      return true;
    }
  }
  return false;
}

::testing::Matcher<assembly> GenotypeIs(std::string expected_gt);

}  // namespace variants

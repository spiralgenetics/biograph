#pragma once

#include "modules/bio_base/biograph.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace variants {

void PrintTo(const assembly& as, std::ostream* os);

class assemble_test : public ::testing::Test {
 protected:
  assemble_test() {
    m_options.bidir_validate_trace_state = 1000;
    m_options.bidir_max_pop_seqset_portion = 1;
    m_options.scaffold = &m_scaffold;
    m_options.report_half_aligned_func = [&](const half_aligned_assembly& ha) {
      std::cout << "Got half-aligned: " << ha << "\n";
    };
    assemble_pipeline_interface::global_set_verify_order(true);
    g_trace_all_assemblies = false;
    reset_assembly_trace();
  }

  pipeline_step_t test_output() {
    return make_unique<assemble_lambda_output>(
        [this](assembly_ptr a) {
          std::cout << "An asssembly was received: " << dump_assembly_and_vars(*a) << "\n";
          if (!(a->right_offset && a->left_offset)) {
            // half-anchored.
            EXPECT_TRUE(a->right_offset || a->left_offset);
            EXPECT_FALSE(a->matches_reference);
            EXPECT_GT(a->seq.size(), 0);
          } else {
            if (a->seq.size() == 0) {
              EXPECT_GT(a->right_offset, a->left_offset) << *a;
            } else {
              EXPECT_GE(a->right_offset, a->left_offset) << *a;
            }
            if (a->right_offset == a->left_offset) {
              EXPECT_GT(a->seq.size(), 0);
            }
          }
          if (a->matches_reference) {
            if (!m_scaffold.empty()) {
              scaffold sub =
                  m_scaffold.subscaffold(a->left_offset, a->right_offset - a->left_offset);
              EXPECT_TRUE(sub.is_simple());
              if (sub.is_simple()) {
                EXPECT_EQ(sub.get_simple(), dna_slice(a->seq));
              }
            }
            // Reference assemblies shouldn't have separate anchors; they're
            // entirely anchors.
            EXPECT_EQ(0, a->left_anchor_len) << *a;
            EXPECT_EQ(0, a->right_anchor_len) << *a;
            EXPECT_EQ(a->seq.size(), a->right_offset - a->left_offset) << *a;
          }
          if (a->matches_reference) {
            m_ref_assemblies.emplace_back(*a);
          } else {
            m_non_ref_assemblies.emplace_back(*a);
          }
          m_assemblies.emplace_back(*a);
          for (const auto& tag : a->tags) {
            m_tag_assemblies[tag].emplace_back(*a);
          }
        },
        "test_output");
  }

  void rev_asm(assembly& a) {
    std::swap(a.left_offset, a.right_offset);
    a.left_offset = m_scaffold.end_pos() - a.left_offset;
    a.right_offset = m_scaffold.end_pos() - a.right_offset;
    a.seq = a.seq.rev_comp();

    if (a.edge_coverage) {
      std::swap(a.edge_coverage->variant_start, a.edge_coverage->variant_end);
      std::swap(a.edge_coverage->reference_start, a.edge_coverage->reference_end);

      auto& ec = *a.edge_coverage;

      for (auto* collection : {&ec.variant_start, &ec.variant_end, &ec.interior,
                               &ec.reference_start, &ec.reference_end}) {
        read_id_set reversed;
        for (uint32_t read_id : *collection) {
          reversed.insert(m_options.readmap->get_rev_comp(read_id));
        }
        std::swap(*collection, reversed);
      }
    }
  }

  void reverse_found_assemblies() {
    for (auto* collection : {&m_assemblies, &m_ref_assemblies, &m_non_ref_assemblies}) {
      for (auto& a : *collection) {
        rev_asm(a);
      }
    }
    for (auto& tagged : m_tag_assemblies) {
      for (auto& a : tagged.second) {
        rev_asm(a);
      }
    }
  }

  void use_biograph(const std::string& bg_dir) {
    m_biograph.emplace(bg_dir);
    m_seqset = m_biograph->get_seqset();
    m_readmap = m_biograph->open_readmap();
    m_options.seqset = m_seqset.get();
    m_options.readmap = m_readmap.get();
  }
  void use_reference(const std::string& reference_dir, const std::string& scaffold_name) {
    m_reference.emplace("", reference_dir);
    m_options.ref = &m_reference.get();
    m_scaffold = trace_ref::ref_to_scaffold(m_options.ref, scaffold_name);
    m_options.scaffold = &m_scaffold;
  }

  void use_ref_parts(std::vector<std::pair<aoffset_t, dna_sequence>> parts) {
    CHECK(m_scaffold.empty());
    aoffset_t last_pos = 0;
    for (const auto& i : parts) {
      EXPECT_GE(i.first, last_pos);
      m_scaffold.add(i.first, i.second);
      last_pos = i.first + i.second.size();
    }

    for (const auto& i : parts) {
      CHECK_EQ(i.second, get_ref_part_seq(i.first, i.second.size())) << " at " << i.first;
    }
  }

  dna_sequence get_ref_part_seq(aoffset_t offset, int len) {
    std::cout << "Getting part seq for " << offset << ", len " << len << "\n";
    scaffold sub = m_scaffold.subscaffold(offset, len);
    CHECK(sub.is_simple());
    dna_slice slice = sub.get_simple();
    return dna_sequence(slice.begin(), slice.end());
  }

  void expect_sorted(const assembly::ordering_t& sort_order) {
    // Verify sorted order
    for (unsigned i = 1; i < m_assemblies.size(); ++i) {
      EXPECT_FALSE(sort_order(m_assemblies[i], m_assemblies[i - 1]))
          << m_assemblies[i - 1] << " should not be before " << m_assemblies[i];
    }
  }

  void use_reads(const std::vector<dna_sequence>& reads) { use_paired_reads({}, reads); }

  uint32_t read_id_for_seq(dna_sequence seq) {
    seqset_range r = m_seqset->find(seq);
    CHECK(r.valid()) << seq;
    CHECK(r.is_seqset_entry()) << seq;
    auto reads = m_readmap->entry_to_index(r.begin());
    for (uint32_t read_id = reads.first; read_id != reads.second; ++read_id) {
      if (m_readmap->get_readlength(read_id) == int(seq.size())) {
        return read_id;
      }
    }
    LOG(FATAL) << "Unable to find read id for sequence: " << seq;
    return std::numeric_limits<uint32_t>::max();
  }

  void use_paired_reads(const std::vector<std::pair<dna_sequence, dna_sequence>>& paired_reads,
                        std::vector<dna_sequence> unpaired_reads = {}) {
    std::vector<dna_sequence> all_reads = unpaired_reads;
    for (const auto& p : paired_reads) {
      all_reads.push_back(p.first);
      all_reads.push_back(p.second);
    }

    m_seqset = seqset_for_reads(all_reads);
    m_readmap = readmap_for_reads(m_seqset, paired_reads, unpaired_reads);

    m_options.seqset = m_seqset.get();
    m_options.readmap = m_readmap.get();
  }

  assemble_options m_options;
  scaffold m_scaffold;
  std::vector<assembly> m_assemblies, m_ref_assemblies, m_non_ref_assemblies;
  std::map<std::string /* tag */, std::vector<assembly>> m_tag_assemblies;
  std::shared_ptr<seqset> m_seqset;
  std::shared_ptr<readmap> m_readmap;
  boost::optional<biograph> m_biograph;
  boost::optional<reference> m_reference;
};

MATCHER_P3(AssemblyIs, left_offset, expected_seq, right_offset, "") {
  return !arg.matches_reference &&                                              //
         bool(arg.left_offset) == bool(optional_aoffset(left_offset)) &&        //
         bool(arg.right_offset) == bool(optional_aoffset(right_offset)) &&      //
         arg.left_offset == left_offset && arg.right_offset == right_offset &&  //
         arg.seq == expected_seq;
}

MATCHER_P2(RefAssemblyIs, left_offset, right_offset, "") {
  return arg.matches_reference && arg.left_offset == aoffset_t(left_offset) &&
         arg.right_offset == aoffset_t(right_offset) && arg.left_anchor_len == 0 &&
         arg.right_anchor_len == 0;
}

class test_pipeline : public pipeline_interface {
 public:
  test_pipeline(pipeline_interface* p) : m_pipeline(p) {}

  pipeline_step_t make_parallel_input() override { return m_pipeline->make_parallel_input(); }

 private:
  pipeline_interface* m_pipeline = nullptr;
};

class test_scaffold_pipeline : public scaffold_pipeline_interface {
 public:
  test_scaffold_pipeline(const std::string& scaffold_name, pipeline_interface* p)
      : m_scaffold_name(scaffold_name), m_pipeline(p) {}

  std::unique_ptr<pipeline_interface> pipeline_for_scaffold(
      const assemble_options& options, const std::string& scaffold_name) override {
    EXPECT_EQ(scaffold_name, m_scaffold_name);
    return make_unique<test_pipeline>(m_pipeline);
  }

 private:
  std::string m_scaffold_name;
  pipeline_interface* m_pipeline = nullptr;
};

std::vector<dna_sequence> reads_for_seq(dna_sequence seq, unsigned read_length,
                                        unsigned read_distance);

namespace coverage_testutil {

class coverage_constructor {
 public:
  operator std::vector<int>() const { return m_cov; }

  coverage_constructor operator+(const coverage_constructor& rhs) const {
    coverage_constructor result = *this;
    result.m_cov.insert(result.m_cov.end(), rhs.m_cov.begin(), rhs.m_cov.end());
    return result;
  }

  coverage_constructor operator+(int n) const {
    coverage_constructor result = *this;
    result.m_cov.push_back(n);
    return result;
  }
  friend coverage_constructor operator+(int lhs, const coverage_constructor& rhs) {
    return coverage_constructor(lhs) + rhs;
  }

  coverage_constructor operator+(const std::vector<int>& rhs) const {
    coverage_constructor result = *this;
    result.m_cov.insert(result.m_cov.end(), rhs.begin(), rhs.end());
    return result;
  }
  friend coverage_constructor operator+(const std::vector<int>& lhs,
                                        const coverage_constructor& rhs) {
    return coverage_constructor(lhs) + rhs;
  }

  coverage_constructor(const std::vector<int>& cov) : m_cov(cov) {}
  coverage_constructor(int cov) : m_cov({cov}) {}

 private:
  std::vector<int> m_cov;
};

MATCHER_P2(CoverageIsInternal, cov, cov_type_arg, "") {
  std::string cov_type_str = cov_type_arg;
  CHECK(cov_type_str == "coverage" || cov_type_str == "pair_coverage");

  const std::vector<int>& actual = cov_type_str == "coverage" ? arg.coverage : arg.pair_coverage;
  if (cov == actual) {
    *result_listener << cov_type_str << " matches\n";
    return true;
  }

  bool need_spacing = false;
  for (int c : cov) {
    if (c > 10) {
      need_spacing = true;
    }
  }
  for (int c : actual) {
    if (c > 10) {
      need_spacing = true;
    }
  }

  *result_listener << "\n" << cov_type_str << " expected:\n";
  int pos = 0;
  for (int c : cov) {
    if ((pos % 10) <= 1) {
      *result_listener << ".";
    } else if (need_spacing) {
      *result_listener << " ";
    }
    ++pos;
    *result_listener << c;
  }
  *result_listener << "\nActual:\n";
  pos = 0;
  for (int c : actual) {
    if ((pos % 10) <= 1) {
      *result_listener << ".";
    } else if (need_spacing) {
      *result_listener << " ";
    }
    ++pos;
    *result_listener << c;
  }
  *result_listener << "\n";

  return false;
}

testing::Matcher<assembly> CoverageIs(const coverage_constructor& cov) {
  return CoverageIsInternal(std::vector<int>(cov), "coverage");
  //, Field(&assembly::coverage, ElementsAreArray(std::vector<int>(cov))));
}

testing::Matcher<assembly> PairCoverageIs(const coverage_constructor& cov) {
  return CoverageIsInternal(std::vector<int>(cov), "pair_coverage");
  //, Field(&assembly::coverage, ElementsAreArray(std::vector<int>(cov))));
}

coverage_constructor rpt(int qty, int n) {
  std::vector<int> result;
  result.resize(qty, n);
  return result;
}

coverage_constructor over(std::string tseq_str, int n) {
  return rpt(dna_testutil::tseq(tseq_str).size() - 1, n);
}

coverage_constructor cov(int n) { return rpt(1, n); }

}  // namespace coverage_testutil

testing::Matcher<const assembly&> AssemblyIdIs(const testing::Matcher<size_t>& id_matcher) {
  return testing::AllOf(testing::Field(&assembly::assembly_id, id_matcher),
                        testing::Field(&assembly::merged_assembly_ids, testing::IsEmpty()));
}

}  // namespace variants

#include "modules/variants/big_assemble_testutil.h"

#include "modules/io/parallel.h"
#include "modules/variants/align.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/calc_coverage.h"
#include "modules/variants/dedup.h"
#include "modules/variants/discovery/state.h"
#include "modules/variants/genotype.h"
#include "modules/variants/normalize.h"
#include "modules/variants/ploid_limit.h"
#include "modules/variants/simple_genotype_filter.h"
#include "modules/variants/trace_ref.h"

#include <fstream>

namespace variants {

using namespace ::testing;

std::vector<std::string> big_assemble_test::g_search_path;
boost::optional<reference> big_assemble_test::m_ref;
std::shared_ptr<seqset> big_assemble_test::m_seqset;
boost::optional<readmap> big_assemble_test::m_readmap;
std::string big_assemble_test::m_refmap_path;
boost::optional<ref_map> big_assemble_test::m_rmap;
std::string big_assemble_test::m_cur_biograph_dir;

void big_assemble_test::init_search_path() {
  CHECK(g_search_path.empty());
  if (getenv("HOME")) {
    g_search_path.push_back(std::string(getenv("HOME")) + std::string("/datasets"));
  }
  if (getenv("USER")) {
    g_search_path.push_back("/home/" + std::string(getenv("USER")) + "/datasets");
  }
  for (const auto& dir : {"/scratch", "/share/datasets/tinyhuman-rand",
                          "/share/datasets/HG002", "/share/datasets/HG001"}) {
    g_search_path.push_back(dir);
  }
}

void big_assemble_test::run_vcf_test(const std::string& scaffold_name,
                                     const std::string& vcf_start_offset /* 1-based offset */,
                                     dna_sequence_matcher ref_bases, dna_sequence_matcher alt_1,
                                     const std::string& gt_1) {
  run_vcf_test_internal(scaffold_name, vcf_start_offset, ref_bases, alt_1, gt_1, boost::none,
                        boost::none);
}
void big_assemble_test::run_vcf_test(const std::string& scaffold_name,
                                     const std::string& vcf_start_offset /* 1-based offset */,
                                     dna_sequence_matcher ref_bases, dna_sequence_matcher alt_1,
                                     const std::string& gt_1, dna_sequence_matcher alt_2,
                                     const std::string& gt_2) {
  run_vcf_test_internal(scaffold_name, vcf_start_offset, ref_bases, alt_1, gt_1, alt_2, gt_2);
}

void big_assemble_test::add_trace(aoffset_t left_offset, aoffset_t right_offset, dna_sequence seq) {
  m_traces.emplace(left_offset, right_offset, seq);
}

void big_assemble_test::run_vcf_test_internal(
    const std::string& scaffold_name, const std::string& vcf_start_offset /* 1-based offset */,
    dna_sequence_matcher ref_bases, dna_sequence_matcher alt_1, const std::string& gt_1,
    boost::optional<dna_sequence_matcher> alt_2, boost::optional<std::string> gt_2) {
  aoffset_t left_offset = std::stol(vcf_start_offset) - 1;
  aoffset_t right_offset = left_offset + ref_bases.size();

  if (m_trace_enabled) {
    m_options.report_bidir_initialized_func = [this](discovery::state* st) {
      for (const auto& t : m_traces) {
        aoffset_t left_offset = std::get<0>(t);
        aoffset_t right_offset = std::get<1>(t);
        CHECK_GT(right_offset, left_offset);
        const auto& seq = std::get<2>(t);
        std::cout << "Adding bidir trace for sequence '" << seq << "' left_offset=" << left_offset
                  << "right_offset= " << right_offset << "\n"
                  << "Ref before: " << m_scaffold.subscaffold_str(left_offset - 200, 200) << "\n"
                  << "Ref during: "
                  << m_scaffold.subscaffold_str(left_offset, right_offset - left_offset) << "\n"
                  << "Ref after: " << m_scaffold.subscaffold_str(right_offset, 200) << "\n";

        st->add_trace_for_variant(left_offset, right_offset, seq);

        constexpr bool k_show_coverage = false;
        if (k_show_coverage) {
          auto ref_around_scaf = m_scaffold.subscaffold(
              left_offset - m_options.seqset->max_read_len(),
              (right_offset - left_offset) + 2 * m_options.seqset->max_read_len());
          if (ref_around_scaf.is_simple()) {
            auto cov = m_options.readmap->approx_coverage(ref_around_scaf.get_simple());
            std::cout << "Ref coverage: " << dump_coverage(cov) << "\n";
          }

          auto left_ref_scaf = m_scaffold.subscaffold(
              left_offset - m_options.seqset->max_read_len(), m_options.seqset->max_read_len());
          auto right_ref_scaf =
              m_scaffold.subscaffold(right_offset, m_options.seqset->max_read_len());
          if (left_ref_scaf.is_simple() && right_ref_scaf.is_simple()) {
            dna_sequence var_cov_seq;
            var_cov_seq += left_ref_scaf.get_simple();
            var_cov_seq += seq;
            var_cov_seq += right_ref_scaf.get_simple();
            auto cov = m_options.readmap->approx_coverage(var_cov_seq);
            std::cout << "Var coverage: " << dump_coverage(cov) << "\n";
          }
        }
      }
    };
  }

  ASSERT_GE(left_offset, 0);

  aoffset_t around_len = std::max<aoffset_t>(100, m_options.readmap->max_read_len());
  aoffset_t around_left = left_offset - around_len;
  if (around_left < 0) {
    around_left = 0;
  }

  aoffset_t around_right = right_offset + around_len;

  std::cout << "Running vcf test for [" << left_offset << ", " << right_offset << "):\n";
  std::cout << "Alt 1, " << gt_1 << ": " << alt_1 << "\n";

  ASSERT_THAT(gt_1, AnyOf(Eq("0/1"), Eq("1/1"), Eq("0/0")))
      << "Genotype must be 0/1 or 1/1 for humans (or 0/0 to ignore genotyping)";
  if (alt_2) {
    ASSERT_THAT(gt_1, Eq("0/1")) << "Genotype must be 0/1 if compound heterozygous";
    ASSERT_THAT(*gt_2, Eq("0/1")) << "Genotype must be 0/1 if compound heterozygous";
  }

  if (alt_2) {
    ASSERT_TRUE(gt_2) << "2nd alt present, but missing 2nd genotype";
  } else {
    ASSERT_FALSE(gt_2) << "No 2nd alt, but 2nd genotype?";
  }

  if (alt_2) {
    std::cout << "Alt 2, " << *gt_2 << ": " << *alt_2 << "\n";
  }

  m_interesting_left_offset = left_offset;
  m_interesting_right_offset = right_offset;
  m_call_pos = (left_offset + right_offset) / 2;
  select_scaffold(scaffold_name);
  if (m_trace_enabled) {
    for (aoffset_t i = left_offset - m_call_around_len; i != left_offset + m_call_around_len; ++i) {
      add_offset_trace(i);
    }
    for (aoffset_t i = right_offset - m_call_around_len; i != right_offset + m_call_around_len;
         ++i) {
      add_offset_trace(m_options.scaffold->end_pos() - i);
    }

    // bidir tracing:
    auto simple1 = alt_1.get_simple();
    if (simple1) {
      add_trace(left_offset, right_offset, *simple1);
    }
    if (alt_2) {
      auto simple2 = alt_2->get_simple();
      if (simple2) {
        add_trace(left_offset, right_offset, *simple2);
      }
    }
  }
  call_region_internal(scaffold_name, left_offset - m_call_around_len,
                       right_offset + m_call_around_len);

  std::string ref_before = m_scaffold.subscaffold_str(around_left, left_offset - around_left);
  std::string ref_after = m_scaffold.subscaffold_str(right_offset, around_right - right_offset);
  std::cout << "Reference region before:\n";
  std::cout << ref_before << "\n";
  std::cout << "Reference region during:";
  std::cout << m_scaffold.subscaffold_str(left_offset, right_offset - left_offset) << "\n";
  std::cout << "Reference region after:\n";
  std::cout << ref_after << "\n";
  std::cout << "\n";

  std::cout << "Expected alt1 sequence:\n"
            << ref_before << " " << alt_1 << " " << ref_after << "\n";
  if (alt_2) {
    std::cout << "Expected alt2 sequence:\n"
              << ref_before << " " << *alt_2 << " " << ref_after << "\n";
  }

  ASSERT_THAT(m_scaffold.subscaffold_str(left_offset, right_offset - left_offset),
              ref_bases.matcher())
      << "Reference does not match what test thinks it is";

  constexpr bool k_check_genotype = false;
  if (k_check_genotype) {
    if (alt_2) {
      EXPECT_THAT(m_assemblies,
                  Contains(AllOf(VariantAt(left_offset, right_offset - left_offset, alt_1),
                                 GenotypeIs(gt_1))))
          << "Alt-1 not called";
      EXPECT_THAT(m_assemblies,
                  Contains(AllOf(VariantAt(left_offset, right_offset - left_offset, *alt_2),
                                 GenotypeIs(*gt_2))))
          << "Alt-2 not called";
    } else {
      EXPECT_THAT(m_assemblies,
                  Contains(AllOf(VariantAt(left_offset, right_offset - left_offset, alt_1),
                                 GenotypeIs(gt_1))))
          << "Alt not called";
    }
  } else {
    EXPECT_THAT(m_assemblies, Contains(VariantAt(left_offset, right_offset - left_offset, alt_1)));
    if (alt_2) {
      EXPECT_THAT(m_assemblies,
                  Contains(VariantAt(left_offset, right_offset - left_offset, alt_1)));
    }
  }

  if (HasNonfatalFailure()) {
    // If we have a fail with variants larger than this size, show
    // which large variants we did find.
    constexpr size_t k_sv_threshold = 50;

    if ((right_offset - left_offset) >= aoffset_t(k_sv_threshold) ||
        alt_1.size() >= k_sv_threshold || (alt_2 && alt_2->size() >= k_sv_threshold)) {
      std::cout << dump_sv_assemblies(k_sv_threshold);
    }
  }
}

pipeline_step_t big_assemble_test::test_output() {
  return make_unique<assemble_lambda_output>(
      [this](assembly_ptr a) {
        std::lock_guard<std::mutex> l(m_mu);
        if (a->left_offset <= m_interesting_right_offset &&
            a->right_offset >= m_interesting_left_offset) {
          std::cout << "Detected call: " << dump_assembly_and_vars(*a) << "\n";
          m_assemblies.emplace_back(*a);
          if (a->matches_reference) {
            m_ref_assemblies.emplace_back(*a);
          } else {
            m_var_assemblies.emplace_back(*a);
          }
        }
        std::cout.flush();
      },
      "test_output");
}

void big_assemble_test::use_biograph(std::string bg_dir) {
  if (m_cur_biograph_dir != bg_dir) {
    open_biograph(bg_dir);
    CHECK_EQ(m_cur_biograph_dir, bg_dir);
  }

  m_options.min_overlap = 80;
  m_options.bidir_max_pop_seqset_portion = 1;
  m_options.bidir_validate_trace_state = 1;
  m_options.debug_paths = [](const std::string& dot_contents) {
    std::string filename = "/tmp/path-debug.dot";
    static std::atomic<size_t> next_debug{0};
    filename += ".";
    filename += std::to_string(next_debug.fetch_add(1));
    std::cout << "Writing path debug to " << filename << "\n";
    std::ofstream path_debug(filename);
    path_debug << dot_contents;
  };
  m_options.seqset = m_seqset.get();
  m_options.readmap = &m_readmap.get();
  m_options.ref = &m_ref.get();
  m_options.rmap = &m_rmap.get();
  m_options.report_aligned_assemblies_func = [this](const assemble_options&, const assembly& a) {
    std::cout << "Got aligned assembly: " << a << ", " << a.aligned_variants.size()
              << " variants:\n";
    for (const auto& v : a.aligned_variants) {
      std::cout << "  " << v << "\n";
    }
    m_aligned.push_back(a);
    m_aligned_dot->add_assembly(a);
  };
  m_options.report_half_aligned_func = [](const half_aligned_assembly& ha) {
    std::cout << "Got half-aligned: " << ha << "\n";
  };
  m_options.report_genotype_discard_func = [](const assemble_options&, const assembly& a,
                                              const std::vector<const assembly*>& better) {
    std::cout << "Genotype discarded assembly: " << dump_assembly_and_vars(a) << "\n";
    constexpr bool k_verbose_discard = false;
    if (k_verbose_discard) {
      std::cout << "because of:\n";
      for (const auto& b : better) {
        if (b->matches_reference) {
          std::cout << "BETTER: REFERENCE id=" << b->assembly_id << "\n";
        } else {
          std::cout << "BETTER:" << dump_assembly_and_vars(*b) << "\n";
        }
      }
    }
  };
}

void big_assemble_test::call_at(const std::string& scaffold_name, std::string one_based_pos,
                                aoffset_t read_around_before, aoffset_t read_around_after) {
  std::cout << "Calling around " << one_based_pos << "\n";
  aoffset_t pos = std::stol(one_based_pos);
  CHECK_GT(pos, 0);
  --pos;
  m_call_pos = pos;

  aoffset_t left = pos - read_around_before;
  aoffset_t right = pos + read_around_after;
  if (left < 0) {
    left = 0;
  }
  scaffold s = trace_ref::ref_to_scaffold(&m_ref.get(), scaffold_name);
  std::cout << "Ref before: " << s.subscaffold_str(left, pos - left) << "\n"
            << "    after:  " << s.subscaffold_str(pos, right - pos) << "\n";

  m_flat_call_pos = m_ref->flatten(scaffold_name, pos);
  m_call_ref_it = m_ref->get_dna(m_flat_call_pos);

  m_interesting_left_offset = pos - 15;
  m_interesting_right_offset = pos + 15;

  select_scaffold(scaffold_name);
  call_region_internal(scaffold_name, left, right);
}

void big_assemble_test::set_thorough_trace_options() {
  m_options.max_ambiguous_search_steps = 500;
  m_options.max_search_steps_per_read = 4;
  m_options.trace_ambiguous_ref = true;
  m_options.max_ploids = 50;
  m_options.max_branches_between_pairs = 5;
  m_options.max_ambiguous_bases = 500;
  m_options.max_cost *= 10;
}

pipeline_step_t big_assemble_test::make_parallel_input() {
  return make_unique<assemble_lambda_copy>(
      [this](const assembly& a) {
        if (a.right_offset < m_interesting_left_offset) {
          return;
        }
        if (a.left_offset > m_interesting_right_offset) {
          return;
        }
        std::lock_guard<std::mutex> l(m_mu);
        std::cout << "Detected raw assembly: ";
        PrintTo(a, &std::cout);
        std::cout << a << "\n";
        std::cout << "Corresponding reference: ";
        if (a.right_offset - a.left_offset <= 1000) {
          std::cout << m_scaffold.subscaffold(a.left_offset, a.right_offset - a.left_offset);
        } else {
          std::cout << m_scaffold.subscaffold(a.left_offset, 500) << "..."
                    << m_scaffold.subscaffold(a.right_offset - 500, 500);
        }
        std::cout << "\n";

        constexpr bool k_dump_pair_reads = false;
        if (k_dump_pair_reads) {
          for (auto matches : {&assembly::left_pair_matches, &assembly::right_pair_matches}) {
            for (uint32_t read_id : a.*matches) {
              std::cout << "Read matched pair: "
                        << m_seqset->ctx_entry(m_readmap->index_to_entry(read_id)).sequence()
                        << ", length " << m_readmap->get_readlength(read_id);
              if (m_readmap->has_mate(read_id)) {
                std::cout << "\nExpected pair: "
                          << m_seqset
                                 ->ctx_entry(
                                     m_readmap->index_to_entry(m_readmap->get_mate(read_id)))
                                 .sequence()
                          << "\n";
              } else {
                std::cout << ", unpaired\n";
              }
            }
          }
        }

        m_assembly_dot->add_assembly(a);
      },
      m_pipeline->make_parallel_input(), "big_assemble_raw");
}

void big_assemble_test::call_region(const std::string& scaffold_name, aoffset_t left_offset,
                                    aoffset_t right_offset) {
  m_interesting_left_offset = left_offset;
  m_interesting_right_offset = right_offset;
  aoffset_t around_len = std::max<aoffset_t>(100, m_options.readmap->max_read_len());
  aoffset_t around_left = std::max(0, left_offset - around_len);
  aoffset_t around_right = right_offset + around_len;
  m_call_pos = (left_offset + right_offset) / 2;
  select_scaffold(scaffold_name);
  call_region_internal(scaffold_name, left_offset, right_offset);
  std::cout << "Reference region before:\n";
  std::cout << m_scaffold.subscaffold_str(around_left, left_offset - around_left) << "\n";
  std::cout << "Reference region during:";
  std::cout << m_scaffold.subscaffold_str(left_offset, right_offset - left_offset) << "\n";
  std::cout << "Reference region after:\n";
  std::cout << m_scaffold.subscaffold_str(right_offset, around_right - right_offset) << "\n";
  std::cout << "\n";
}

void big_assemble_test::select_scaffold(const std::string& scaffold_name) {
  m_scaffold = trace_ref::ref_to_scaffold(m_options.ref, scaffold_name);
  m_options.scaffold = &m_scaffold;
}

void big_assemble_test::call_region_internal(const std::string& scaffold_name, aoffset_t start,
                                             aoffset_t limit) {
  m_assembly_dot.emplace(m_scaffold);
  m_aligned_dot.emplace(m_scaffold);
  std::cout << "Scaffold end pos: " << m_scaffold.end_pos() << "\n";
  std::cout << "Call distance from end: " << m_scaffold.end_pos() - m_call_pos << "\n";
  m_pipeline = make_unique<assemble_pipeline>(m_options, test_output());
  m_pipeline->add_standard_variants_pipeline();
  m_options.scaffold = nullptr;
  {
    test_scaffold_pipeline sp(scaffold_name, this);
    {
      trace_ref t(m_options, &sp);
      aoffset_t trace_start = start;
      aoffset_t trace_limit = limit;
      if (m_options.min_pair_depth > 0 || m_options.min_avg_pair_depth > 0) {
        // Make sure we get pair info
        trace_start -= m_options.max_pair_distance;
        trace_limit += m_options.max_pair_distance;
        std::cout << "Adjusting call range to [" << trace_start << ", " << trace_limit
                  << ") to make sure we get pairing data\n";
      }
      t.add_scaffold_range(scaffold_name, trace_start, trace_limit);
      auto old_thread_count = get_thread_count();
      set_thread_count("1");
      m_stats = t.assemble();
      set_thread_count(std::to_string(old_thread_count));
      std::cout << "Assemble stats: " << m_stats << "\n";
    }
    std::cout << "Finishing assembly pipeline\n";
    m_pipeline.reset();
  }
  static size_t g_dot_index = 0;
  {
    std::string assemblies_dot = "/tmp/assemblies.dot." + std::to_string(g_dot_index);
    std::cout << "Writing assembly dot to " << assemblies_dot << "\n";
    std::ofstream f(assemblies_dot);
    f << m_assembly_dot->str();
  }
  {
    std::string aligned_dot = "/tmp/aligned.dot." + std::to_string(g_dot_index);
    std::cout << "Writing aligned dot to " << aligned_dot << "\n";
    std::ofstream f(aligned_dot);
    f << m_aligned_dot->str();
  }
  ++g_dot_index;
}

void big_assemble_test::open_biograph(std::string bg_name) {
  m_rmap.reset();
  m_ref.reset();
  m_readmap.reset();
  m_seqset.reset();
  m_readmap.reset();

  if (g_search_path.empty()) {
    init_search_path();
  }
  std::cout << "BioGraph search path: " << PrintToString(g_search_path) << "\n";
  std::string full_bg_path;
  for (const auto& d : g_search_path) {
    if (boost::filesystem::exists(d + "/" + bg_name)) {
      full_bg_path = d + "/" + bg_name;
      if (boost::filesystem::exists(full_bg_path + ".refmap")) {
        m_refmap_path = full_bg_path + ".refmap";
      }
      break;
    }
  }
  CHECK(!full_bg_path.empty()) << "Could not find " << bg_name << " in "
                               << PrintToString(g_search_path);
  std::cout << "Using " << bg_name << " in " << full_bg_path << "\n";

  m_ref.emplace("", "/reference/hs37d5");
  biograph_dir bgdir(full_bg_path, open_mode::READ_BGDIR);

  std::string seqset_path = bgdir.seqset();
  std::string readmap_path = bgdir.find_readmap("");

  std::cout << "Using seqset: " << seqset_path << " readmap: " << readmap_path << "\n";
  CHECK(!seqset_path.empty());
  CHECK(!readmap_path.empty());
  m_seqset = std::make_shared<seqset>(seqset_path);
  m_readmap.emplace(m_seqset, readmap_path);
  if (m_refmap_path.empty()) {
    m_refmap_path = "/scratch/" + bg_name + ".refmap";
    if (!boost::filesystem::exists(m_refmap_path)) {
      std::cout << "Building refmap\n\n";
      std::cout.flush();
      std::string new_refmap = m_refmap_path + ".new";
      unlink(new_refmap.c_str());
      try {
        {
          spiral_file_create_mmap c(new_refmap);
          m_rmap.emplace(m_seqset.get(), &m_ref.get(), c.create());
          std::cout << "Starting build..\n";
          std::cout.flush();
          m_rmap->build(noisy_progress_handler());
        }
        std::cout << "Refmap build complete\n";
        rename(new_refmap.c_str(), m_refmap_path.c_str());
      } catch (const io_exception& e) {
        std::cout << "Can't write to refmap " << m_refmap_path << ": " << e.what() << "\n";
        std::cout << "Building refmap in memory.\n";
        std::cout.flush();
        m_refmap_path.clear();
        m_rmap.emplace(m_seqset.get(), &m_ref.get());
        m_rmap->build(noisy_progress_handler());
      }
    }
  }
  if (!m_rmap) {
    CHECK(!m_refmap_path.empty());
    std::cout << "Opening refmap: " << m_refmap_path << "\n";
    std::cout.flush();
    spiral_file_open_mmap o(m_refmap_path);
    m_rmap.emplace(m_seqset.get(), &m_ref.get(), o.open());
  }

  m_cur_biograph_dir = bg_name;
}

std::string big_assemble_test::dump_sv_assemblies(int min_size) {
  std::stringstream out;

  out << "Filtering out assemblies smaller than " << min_size << " bases from "
      << m_assemblies.size() << " assemblies:\n";
  int too_small = 0;
  for (const auto& a : m_assemblies) {
    if (a.right_offset - a.left_offset < min_size && int(a.seq.size()) < min_size) {
      too_small++;
    } else {
      out << dump_assembly_and_vars(a) << "\n";
    }
  }
  out << "... plus " << too_small << " assemblies smaller than " << min_size << " bases\n";

  return out.str();
}

::testing::Matcher<assembly> GenotypeIs(std::string expected_gt) {
  if (expected_gt == "0/1") {
    return Field(&assembly::strand_count, 1);
  } else if (expected_gt == "1/1") {
    return Field(&assembly::strand_count, 2);
  } else if (expected_gt == "0/0") {
    // Don't check genotyping
    return _;
  } else {
    LOG(FATAL) << "Unknown genotype " << expected_gt;
    return Field(&assembly::strand_count, 0);
  }
}

}  // namespace variants

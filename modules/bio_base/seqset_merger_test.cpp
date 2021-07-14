#include "modules/bio_base/seqset_merger.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/make_mergemap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/io/parallel.h"
#include "modules/io/spiral_file.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/test/test_coverage.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;

class seqset_merger_test : public Test {
 public:
  void add_input_seqs(const std::vector<dna_sequence>& seqs);
  void verify();
  void run_coverage_pass(unsigned random_seed);
  void reset();

 private:
  void make_mergemaps();
  void merge();

  std::vector<std::unique_ptr<seqset_file>> m_seqsets;
  std::vector<std::unique_ptr<seqset_flat>> m_flats;
  std::vector<std::unique_ptr<seqset_mergemap>> m_mergemaps;

  std::unique_ptr<spiral_file_create_mem> m_merge_create;
  std::unique_ptr<seqset> m_merged_seqset;
};

void seqset_merger_test::reset() {
  m_merge_create.reset();
  m_merged_seqset.reset();
  m_mergemaps.clear();
  m_flats.clear();
  m_seqsets.clear();
}

void seqset_merger_test::add_input_seqs(const std::vector<dna_sequence>& seqs) {
  m_seqsets.emplace_back(seqset_for_reads(seqs));
  m_flats.emplace_back(
      seqset_flat_for_seqset(&m_seqsets.back().get()->get_seqset()));
}

void seqset_merger_test::make_mergemaps() {
  CHECK(!m_merge_create);
  m_merge_create.reset(new spiral_file_create_mem());

  std::vector<const seqset_flat*> flats;
  for (const auto& flat : m_flats) {
    flats.push_back(flat.get());
  }

  make_mergemap make_mm(flats);
  make_mm.build();

  for (unsigned i = 0; i < m_seqsets.size(); ++i) {
    spiral_file_create_mem c;
    seqset_mergemap_builder mm(c.create(), m_seqsets[i]->get_seqset().uuid(),
                               m_merge_create->uuid(),
                               make_mm.total_merged_entries());
    make_mm.fill_mergemap(i, &mm);
    spiral_file_mem_storage encoded = c.close();
    spiral_file_open_mem o(encoded);
    m_mergemaps.emplace_back(new seqset_mergemap(o.open()));
  }
}

void seqset_merger_test::merge() {
  std::vector<const seqset_flat*> flats;
  for (const auto& flat : m_flats) {
    flats.push_back(flat.get());
  }

  std::vector<const seqset_mergemap*> mergemaps;
  for (const auto& mergemap : m_mergemaps) {
    mergemaps.push_back(mergemap.get());
  }

  seqset_merger merger(flats, mergemaps);
  merger.build(m_merge_create->create());
  spiral_file_mem_storage encoded = m_merge_create->close();

  spiral_file_open_mem o(encoded);
  m_merged_seqset.reset(new seqset(o.open()));
}

void seqset_merger_test::verify() {
  make_mergemaps();
  merge();

  std::set<dna_sequence> expected_seqs;

  for (const auto& ss_f : m_seqsets) {
    const auto& ss = ss_f->get_seqset();
    for (size_t i = 0; i < ss.size(); ++i) {
      expected_seqs.insert(ss.ctx_entry(i).sequence());
    }
  }

  std::set<dna_sequence>::iterator it, next;
  for (it = expected_seqs.begin(); it != expected_seqs.end(); it = next) {
    next = it;
    ++next;

    if (next != expected_seqs.end()) {
      // Remove any prefixes of other entries.
      if (next->size() > it->size() && next->subseq(0, it->size()) == *it) {
        expected_seqs.erase(it);
      }
    }
  }

  std::set<dna_sequence> actual_seqs;
  for (size_t i = 0; i < m_merged_seqset->size(); ++i) {
    actual_seqs.insert(m_merged_seqset->ctx_entry(i).sequence());
  }

  EXPECT_THAT(actual_seqs, ContainerEq(expected_seqs));
}

using namespace dna_testutil;

TEST_F(seqset_merger_test, single_simple) {
  add_input_seqs({tseq("abc"), tseq("de")});
  verify();
}

TEST_F(seqset_merger_test, merge2) {
  add_input_seqs({tseq("abc"), tseq("cde")});
  add_input_seqs({tseq("abc"), tseq("efg")});
  verify();
}

void seqset_merger_test::run_coverage_pass(unsigned rand_seed) {
  std::cerr << "Running coverage pass with random seed: " << rand_seed << "\n";
  std::mt19937 random_source(rand_seed);

  auto orig_parallel_splits = g_parallel_splits;
  std::uniform_int_distribution<unsigned> parallel_splits_gen(1, 100);
  g_parallel_splits = parallel_splits_gen(random_source);

  std::uniform_int_distribution<unsigned> merge_parts_gen(1, 5);
  std::uniform_int_distribution<unsigned> seq_len_gen(5, 20);
  std::uniform_int_distribution<unsigned> seq_count_gen(10, 20);

  unsigned n_parts = merge_parts_gen(random_source);
  std::cerr << "Generating " << n_parts << " merge parts\n";
  for (unsigned i = 0; i < n_parts; ++i) {
    std::vector<dna_sequence> part_seqs;
    size_t num_seqs = seq_count_gen(random_source);
    std::cerr << "Generating " << num_seqs << " random seqs "
              << ":\n";
    for (size_t j = 0; j < num_seqs; ++j) {
      part_seqs.push_back(
          rand_dna_sequence(random_source, seq_len_gen(random_source)));
      std::cerr << " " << part_seqs.back().as_string() << "\n";
    }
    add_input_seqs(part_seqs);
  }
  verify();
  reset();
  g_parallel_splits = orig_parallel_splits;
}

TEST_F(seqset_merger_test, coverage) {
  scoped_test_coverage cov;

  std::mt19937 seed_source(time(0));

  std::string last_missing;
  while (!cov.missing("seqset_merger").empty()) {
    std::string missing =
        ::testing::PrintToString(cov.missing("seqset_merger"));
    if (missing != last_missing) {
      last_missing = missing;
      std::cerr << "Missing coverage: " << missing << "\n";
    }

    run_coverage_pass(seed_source());
  }
}

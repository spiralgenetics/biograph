#include "modules/bio_base/make_mergemap.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/io/uuid.h"
#include "modules/test/test_coverage.h"

using namespace testing;
using namespace dna_testutil;

namespace {

void merge_and_verify(const std::vector<std::vector<dna_sequence>>& parts) {
  unsigned n_merge = parts.size();
  CHECK_GE(n_merge, 1);

  std::vector<std::unique_ptr<seqset_file>> ss;
  std::vector<std::unique_ptr<seqset_flat>> ss_flat;
  std::vector<const seqset_flat*> ss_flat_ptrs;

  for (unsigned i = 0; i < n_merge; ++i) {
    ss.emplace_back(seqset_for_reads(parts[i]));
    ss_flat.emplace_back(seqset_flat_for_seqset(&ss[i]->get_seqset()));
    ss_flat_ptrs.emplace_back(ss_flat[i].get());

    // dump_seqset("P#" + std::to_string(i) + ": ", ss[i]->get_seqset());
  }

  std::vector<dna_sequence> all_reads;
  for (const auto& part : parts) {
    for (const auto& seq : part) {
      all_reads.push_back(seq);
    }
  }
  std::unique_ptr<seqset_file> merged_ss_file = seqset_for_reads(all_reads);
  const seqset& merged_ss = merged_ss_file->get_seqset();
  std::unique_ptr<seqset_flat> merged_ss_flat =
      seqset_flat_for_seqset(&merged_ss);

  // dump_seqset("M: ", merged_ss);

  make_mergemap counter(ss_flat_ptrs);
  counter.build();

  EXPECT_EQ(merged_ss.size(), counter.total_merged_entries());

  std::string merged_uuid = make_uuid();

  std::vector<spiral_file_mem_storage> mergemaps;
  for (unsigned i = 0; i < n_merge; ++i) {
    spiral_file_create_mem c;
    seqset_mergemap_builder mm(c.create(), ss[i]->get_seqset().uuid(),
                               merged_uuid, counter.total_merged_entries());
    counter.fill_mergemap(i, &mm);
    mm.finalize();
    mergemaps.push_back(c.close());
  }

  for (unsigned i = 0; i < n_merge; ++i) {
    spiral_file_open_mem o(mergemaps[i]);
    seqset_mergemap mm(o.open());
    const bitcount& bc = mm.get_bitcount();

    EXPECT_EQ(ss_flat[i]->size(), bc.total_bits());
    EXPECT_EQ(ss[i]->get_seqset().uuid(), mm.metadata().orig_seqset_uuid);
    EXPECT_EQ(merged_uuid, mm.metadata().merged_seqset_uuid);

    std::set<size_t> wrongly_set_bits;
    for (size_t idx = 0; idx < merged_ss.size(); ++idx) {
      if (bc.get(idx)) {
        wrongly_set_bits.insert(idx);
      }
    }

    for (size_t part_idx = 0; part_idx < ss_flat[i]->size(); part_idx++) {
      dna_slice slice = ss_flat[i]->get(part_idx);
      seqset_range entry = merged_ss.find(slice);
      if (!entry.valid()) {
        continue;
      }

      EXPECT_TRUE(bc.get(entry.begin()))
          << "part " << i << ": part entry " << part_idx
          << ": merge entry: " << entry.begin() << "-" << entry.end() << ": "
          << entry.sequence().as_string();
      wrongly_set_bits.erase(entry.begin());
    }
    EXPECT_THAT(wrongly_set_bits, IsEmpty()) << "part " << i;
  }
}

void run_coverage_pass(unsigned rand_seed) {
  std::cerr << "Running coverage pass with random seed: " << rand_seed << "\n";
  std::mt19937 random_source(rand_seed);

  std::vector<std::vector<dna_sequence>> seqs;
  seqs.resize(std::uniform_int_distribution<unsigned>(1, 15)(random_source));
  std::uniform_int_distribution<unsigned> seq_len_gen(5, 20);
  std::uniform_int_distribution<unsigned> seq_count_gen(10, 20);
  for (std::vector<dna_sequence>& part_seqs : seqs) {
    size_t num_seqs = seq_count_gen(random_source);
    std::cerr << "Generating " << num_seqs << " random seqs "
              << ":\n";
    for (size_t i = 0; i < num_seqs; ++i) {
      part_seqs.push_back(
          rand_dna_sequence(random_source, seq_len_gen(random_source)));
      std::cerr << " " << part_seqs.back().as_string() << "\n";
    }
  }
  merge_and_verify(seqs);
}

}  // namespace

TEST(make_mergemap, merge_single) {
  std::vector<dna_sequence> seqs = {tseq("ab"), tseq("bc"), tseq("cd"),
                                    tseq("be")};

  merge_and_verify({seqs});
}

TEST(make_mergemap, merge_two) {
  std::vector<dna_sequence> seqs1 = {tseq("ab"), tseq("bc"), tseq("cd"),
                                     tseq("be")};
  std::vector<dna_sequence> seqs2 = {tseq("AB"), tseq("BC"), tseq("CD"),
                                     tseq("BE")};

  merge_and_verify({seqs1, seqs2});
}

TEST(make_mergemap, coverage) {
  scoped_test_coverage cov;

  std::mt19937 seed_source(time(0));

  std::string last_missing;
  while (!cov.missing("make_mergemap").empty()) {
    std::string missing =
        ::testing::PrintToString(cov.missing("make_mergemap"));
    if (missing != last_missing) {
      last_missing = missing;
      std::cerr << "Missing coverage: " << missing << "\n";
    }

    run_coverage_pass(seed_source());
  }
}

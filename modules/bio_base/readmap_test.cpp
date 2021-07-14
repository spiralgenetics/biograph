#include "gtest/gtest.h"

#include "modules/io/file_io.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/io/json_transfer.h"
#include "modules/mapred/path.h"
#include "modules/mapred/output_stream.h"

#include "modules/bio_base/seqset.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset_mergemap.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/make_mergemap.h"
#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_base/seqset_merger.h"
#include "modules/bio_base/biograph_dir.h"

#include "modules/bio_mapred/make_readmap.h"
#include "modules/bio_mapred/flatten_seqset.h"
#include "modules/bio_mapred/merge_flat_seqset.h"
#include "modules/test/test_utils.h"
#include "modules/bio_base/biograph.h"

#include <gmock/gmock.h>

using namespace testing;
using namespace dna_testutil;

manifest read_corrected_reads(const std::string& manifest_path) {
  std::string serialized_manifest{slurp_file(manifest_path)};
  return inline_json_deserialize<manifest>(serialized_manifest);
}

TEST(readmap, get_bit) {
  biograph bg("datasets/hiv/biograph/ERR732131.bg");

  std::shared_ptr<seqset> the_seqset = bg.get_seqset();
  std::shared_ptr<readmap> the_readmap = bg.open_readmap();

  SPLOG("SEQSET has %lu entries", the_seqset->size());

  for (auto entry_id = 0UL; entry_id < the_seqset->size(); entry_id++) {
    seqset_range entry_range = the_seqset->ctx_entry(entry_id);

    // NOTE: this test only works if autotrimming is disabled during biograph creation.
    ASSERT_EQ(entry_range.is_full_read(*the_readmap) == true, the_readmap->get_bit(entry_id));
  }
}

TEST(readmap, read_props_tests) {
  // Few random type reads to look for
  std::vector<std::vector<dna_sequence>> test_pairs = {
      {tseq("READ1"), tseq("ANOTHER1")},
      {tseq("NEWB"), tseq("BROTHER")},
      {tseq("SOLO"), tseq("")},
      {tseq("PREFIXread"), tseq("PREFIXmate")},
      {tseq("readSUFFIX"), tseq("mateSUFFIX")},
      {tseq("PREFIXreadSUFFIX"), tseq("PREFIXmateSUFFIX")},
      {tseq("read"), tseq("mate")},
      {tseq("XreadS"), tseq("XmateS")},
      //Test case from NIS that would fail..
      //{tseq("1ABCD"), tseq("")}, But I'd need to do the entry_read.end/entry_read.begin
      // {tseq("2ABCD"), tseq("")},
      // {tseq("3ABCD"), tseq("")},
      // {tseq("ABCD"), tseq("")},
  };
  // let's open up a biograph
  auto bgfile = biograph_for_reads(test_pairs);
  std::shared_ptr<seqset> the_seqset = std::move(bgfile.first);
  std::unique_ptr<readmap> the_readmap = std::move(bgfile.second);

  auto props_test = [&](const dna_sequence& read, const dna_sequence& pair,
                        bool fwd = true) {
    // Find the entry
    auto entry_read = the_seqset->find(read);
    ASSERT_TRUE(entry_read.valid());
    ASSERT_EQ(entry_read.end() - entry_read.begin(), 1);

    // Turn the entry into a read index
    ASSERT_TRUE(the_readmap.get()->get_bit(entry_read.begin()));
    auto read_idx_p =
        the_readmap.get()->entry_to_index(entry_read.begin());

    uint32_t read_idx;
    if (read_idx_p.second - read_idx_p.first != 1) {
      // Multiple reads in this entry, so pick the right one
      bool set = false;  // Did we find a read to use
      for (auto i = read_idx_p.first; i < read_idx_p.second; i++) {
        if ((unsigned)the_readmap.get()->get_readlength(i) == read.size()) {
          read_idx = i;
          set = true;
        }
      }
      ASSERT_TRUE(set);
    } else {
      read_idx = read_idx_p.first;
    }
    // Lookup the properties on the read index
    ASSERT_EQ(the_readmap.get()->get_is_forward(read_idx), fwd); //skip real quick
    ASSERT_EQ(the_readmap.get()->get_readlength(read_idx), read.size());

    // Work backwards to the entry
    auto entry_id = the_readmap.get()->index_to_entry(read_idx);
    auto entry_ret = the_seqset->read_ctx_entry(
        *the_readmap.get(), read_idx);
    ASSERT_EQ(entry_id, entry_read.begin());
    ASSERT_EQ(read, entry_ret.sequence());

    // Check finding the reverse complement read.  TODO(nils): Make
    // this test more stringent once we store reverse complements with
    // the readmap; the pair's reverse complement should be paired
    // with our reverse complement, even in ambiguous cases.
    auto rc_read_idx = the_readmap.get()->get_rev_comp(read_idx);
    auto rc_ret = the_seqset->read_ctx_entry(
        *the_readmap.get(), rc_read_idx);
    EXPECT_EQ(entry_ret.sequence().rev_comp(), rc_ret.sequence());

    // Pair work
    if (pair.size() != 0) {
      auto entry_pair = the_seqset->find(pair);
      ASSERT_TRUE(entry_pair.valid());
      ASSERT_EQ(entry_pair.end() - entry_pair.begin(), 1);

      ASSERT_TRUE(the_readmap.get()->get_bit(entry_pair.begin()));
      auto pair_idx_p =
          the_readmap.get()->entry_to_index(entry_pair.begin());
      uint32_t pair_idx;
      if (pair_idx_p.second - pair_idx_p.first != 1) {
        // Multiple reads in this entry, so pick the right one
        bool set = false;
        for (auto i = pair_idx_p.first; i < pair_idx_p.second; i++) {
          if ((unsigned)the_readmap.get()->get_readlength(i) == pair.size()) {
            pair_idx = i;
            set = true;
          }
        }
        ASSERT_TRUE(set);
      } else {
        pair_idx = pair_idx_p.first;
      }

      // Lookup the properties on the read index
      EXPECT_EQ(the_readmap.get()->get_is_forward(pair_idx), fwd);
      EXPECT_EQ(the_readmap.get()->get_readlength(pair_idx), pair.size());

      // Pairs point to each other
      ASSERT_TRUE(the_readmap.get()->has_mate(read_idx));
      ASSERT_TRUE(the_readmap.get()->has_mate(pair_idx));
      EXPECT_EQ(the_readmap.get()->get_mate(read_idx), pair_idx) << "read_idx: " << read_idx;
      EXPECT_EQ(the_readmap.get()->get_mate(pair_idx), read_idx) << "pair_idx: " << pair_idx;

      // Work backwards to the entry
      auto entry_pair_id = the_readmap.get()->index_to_entry(pair_idx);
      auto entry_pair_ret =
          the_seqset->read_ctx_entry(*the_readmap.get(), pair_idx);
      ASSERT_EQ(entry_pair_id, entry_pair.begin());//??
      ASSERT_EQ(pair, entry_pair_ret.sequence());
    }
  };

  for (auto pair : test_pairs) {
    props_test(pair[0], pair[1]);
    props_test(pair[0].rev_comp(), pair[1].rev_comp(), false);
  }
}

TEST(readmap, migrate) {
  std::string merged_file_path{make_path("merged_hiv_seqset")};
  std::string migrated_readmap_path{make_path("merged_hiv_readmap")};

  unlink(merged_file_path.c_str());
  unlink(migrated_readmap_path.c_str());

  SPLOG("Loading biograph...");
  biograph bg("datasets/hiv/biograph/ERR732130.bg");
  std::shared_ptr<seqset> original_seqset = bg.get_seqset();
  std::shared_ptr<readmap> original_readmap = bg.open_readmap();

  SPLOG("Merging HIV seqsets...");
  std::vector<std::string> seqsets {
      "datasets/hiv/biograph/ERR381524.bg/seqset", "datasets/hiv/biograph/ERR732129.bg/seqset",
      "datasets/hiv/biograph/ERR732131.bg/seqset", "datasets/hiv/biograph/ERR732132.bg/seqset",
      "datasets/hiv/biograph/ERR732130.bg/seqset"
  };

  SPLOG("Flattening SEQSET...");
  flatten_seqset flattener(seqsets, 32);
  std::multimap<int, std::shared_ptr<scoped_temp_file>> temp_file_map = flattener();

  SPLOG("Merging SEQSET...");
  merge_flat_seqsets()(merged_file_path, temp_file_map, true, 255 /* max read len */);
  auto merged_seqset = std::make_shared<seqset>(merged_file_path);

  SPLOG("Migrating readmap to merged seqset...");
  make_readmap::migrate(*original_seqset, *original_readmap, *merged_seqset,
                        migrated_readmap_path);

  readmap merged_readmap(merged_seqset, migrated_readmap_path);

  for (auto i = 0ULL; i < original_seqset->size(); i++) {
    if (original_readmap->get_bit(i)) {
      seqset_range original_entry = original_seqset->ctx_entry(i);
      ASSERT_TRUE(original_entry.valid());
      ASSERT_EQ(original_entry.size(), 250);

      const auto original_sequence = original_entry.sequence();
      seqset_range merged_range = merged_seqset->find(original_sequence);
      ASSERT_TRUE(merged_range.valid());
      ASSERT_EQ(merged_range.size(), 250);
      ASSERT_EQ(merged_range.begin(), merged_range.end() - 1);
      ASSERT_EQ(merged_range.sequence(), original_sequence);
      ASSERT_TRUE(merged_readmap.get_bit(merged_range.begin()));
    }
  }

  // Walk merged seqset and find in original and match getbit from readmap.
  for (auto i = 0ULL; i < merged_seqset->size(); i++) {
    if (merged_readmap.get_bit(i)) {
      seqset_range merged_range = merged_seqset->ctx_entry(i);
      ASSERT_TRUE(merged_range.valid());
      ASSERT_EQ(merged_range.begin(), merged_range.end() - 1);
      auto merged_seq = merged_range.sequence();
      ASSERT_EQ(merged_seq.size(), 250);
      seqset_range original_range = original_seqset->find(merged_seq);
      if (original_range.valid()) {
        ASSERT_TRUE(original_readmap->get_bit(original_range.begin()));
      }
    }
  }
}


TEST(readmap, fast_migrate) {
  std::vector<std::unique_ptr<seqset_mergemap>> mergemaps;
  std::vector<const seqset_mergemap*> mergemap_ptrs;

  std::string merged_file_path{make_path("merged_hiv_seqset")};
  std::string migrated_readmap_path{make_path("merged_hiv_readmap")};

  unlink(merged_file_path.c_str());
  unlink(migrated_readmap_path.c_str());

  SPLOG("Loading biograph...");
  biograph bg("datasets/hiv/biograph/ERR732131.bg");
  std::shared_ptr<seqset> original_seqset = bg.get_seqset();
  std::shared_ptr<readmap> original_readmap = bg.open_readmap();

  SPLOG("Merging HIV seqsets...");
  std::vector<std::string> seqsets {
      "datasets/hiv/biograph/ERR381524.bg/seqset", "datasets/hiv/biograph/ERR732129.bg/seqset",
      "datasets/hiv/biograph/ERR732131.bg/seqset", "datasets/hiv/biograph/ERR732132.bg/seqset",
      "datasets/hiv/biograph/ERR732130.bg/seqset"
  };

  std::vector<std::shared_ptr<seqset>> ss_files;
  for (const auto& input : seqsets) {
    ss_files.emplace_back(std::make_shared<seqset>(input));
  }

  std::vector<std::unique_ptr<seqset_flat>> flats;
  std::vector<const seqset_flat*> flat_ptrs;
  for (const auto& ss_f : ss_files) {
    flats.emplace_back(seqset_flat_for_seqset(&ss_f->get_seqset()));
    flat_ptrs.emplace_back(flats.back().get());
  }

  {
    spiral_file_create_mmap merge_create(merged_file_path);

    make_mergemap make_mm(flat_ptrs);
    make_mm.build();

    for (unsigned i = 0; i < ss_files.size(); ++i) {
      spiral_file_create_mem c;
      seqset_mergemap_builder mm(c.create(), ss_files[i]->get_seqset().uuid(),
                                 merge_create.uuid(),
                                 make_mm.total_merged_entries());
      make_mm.fill_mergemap(i, &mm);
      spiral_file_mem_storage encoded = c.close();
      spiral_file_open_mem o(encoded);
      mergemaps.emplace_back(new seqset_mergemap(o.open()));
      mergemap_ptrs.emplace_back(mergemaps.back().get());
    }

    seqset_merger merger(flat_ptrs, mergemap_ptrs);
    merger.build(merge_create.create());

    SPLOG("Creating new readmap...");
    spiral_file_create_mmap new_readmap(migrated_readmap_path);

    SPLOG("Starting fast_migrate...");
    make_readmap::fast_migrate(*original_readmap, *mergemaps[2],
                               new_readmap.create());
  }

  SPLOG("Opening merged seqset...");
  auto merged_seqset = std::make_shared<seqset>(merged_file_path);

  SPLOG("Opening merged readmap...");
  readmap merged_readmap(merged_seqset, migrated_readmap_path);

  for (auto i = 0ULL; i < original_seqset->size(); i++) {
    if (original_readmap->get_bit(i)) {
      seqset_range original_entry = original_seqset->ctx_entry(i);
      ASSERT_TRUE(original_entry.valid());
      ASSERT_EQ(original_entry.size(), 250);

      const auto original_sequence = original_entry.sequence();
      seqset_range merged_range = merged_seqset->find(original_sequence);
      ASSERT_TRUE(merged_range.valid());
      ASSERT_EQ(merged_range.size(), 250);
      ASSERT_EQ(merged_range.begin(), merged_range.end() - 1);
      ASSERT_EQ(merged_range.sequence(), original_sequence);
      ASSERT_TRUE(merged_readmap.get_bit(merged_range.begin()));
    }
  }

  // Walk merged seqset and find in original and match getbit from readmap.
  for (auto i = 0ULL; i < merged_seqset->size(); i++) {
    if (merged_readmap.get_bit(i)) {
      seqset_range merged_range = merged_seqset->ctx_entry(i);
      ASSERT_TRUE(merged_range.valid());
      ASSERT_EQ(merged_range.begin(), merged_range.end() - 1);
      auto merged_seq = merged_range.sequence();
      ASSERT_EQ(merged_seq.size(), 250);
      seqset_range original_range = original_seqset->find(merged_seq);
      if (original_range.valid()) {
        ASSERT_TRUE(original_readmap->get_bit(original_range.begin()));
      }
    }
  }
}

std::vector<dna_sequence> get_prefix_read_seqs(const readmap* rm,
                                               const seqset_range& r) {
  std::vector<dna_sequence> result;
  for (const auto& read : rm->get_prefix_reads(r)) {
    seqset_range read_r = read.get_seqset_entry();
    result.push_back(read_r.sequence());
  }
  return result;
}

TEST(readmap, get_prefix_reads) {
  std::vector<std::vector<dna_sequence>> test_reads =  //
      {{tseq("z")},
       {tseq("abcde")},
       {tseq("abcde") + dna_A},
       {tseq("abcde") + dna_A},
       {tseq("abcde") + dna_C},
       {tseq("abcde") + dna_C + tseq("f")},
       {tseq("abcde") + dna_T}};
  auto bgfile = biograph_for_reads(test_reads);
  std::shared_ptr<seqset> the_seqset = std::move(bgfile.first);
  std::unique_ptr<readmap> the_readmap = std::move(bgfile.second);

  seqset_range r;

  r = the_seqset->find(tseq("abcd"));
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r), IsEmpty());

  r = the_seqset->find(tseq("abcde"));
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r), UnorderedElementsAre(tseq("abcde")));

  r = the_seqset->find(tseq("abcde") + dna_A);
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_A, tseq("abcde") + dna_A));

  r = the_seqset->find(tseq("abcde") + dna_C);
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_C));

  r = the_seqset->find(tseq("abcde") + dna_C + tseq("f"));
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_C,
                                   tseq("abcde") + dna_C + tseq("f")));

  r = the_seqset->find(tseq("abcde") + dna_T);
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_T));
}

TEST(readmap, get_prefix_reads2) {
  // Same as get_prefix_reads, but adds some additional seqset entries in the middle that are
  // longer.
  std::vector<std::vector<dna_sequence>> test_reads =  //
      {{tseq("z")},
       {tseq("abcde")},
       {dna_T + tseq("abcde") + dna_A + dna_A},
       {tseq("abcde") + dna_A},
       {tseq("abcde") + dna_A},
       {dna_C + tseq("abcde") + dna_G},
       {tseq("abcde") + dna_C},
       {tseq("abcde") + dna_C + tseq("f")},
       {tseq("abcde") + dna_T},
       {dna_A + tseq("abcde") + dna_T + dna_T}};

  auto bgfile = biograph_for_reads(test_reads);
  std::shared_ptr<seqset> the_seqset = std::move(bgfile.first);
  std::unique_ptr<readmap> the_readmap = std::move(bgfile.second);

  seqset_range r;

  r = the_seqset->find(tseq("a"));
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r), IsEmpty());

  r = the_seqset->find(tseq("abcd"));
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r), IsEmpty());

  r = the_seqset->find(tseq("abcde"));
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r), UnorderedElementsAre(tseq("abcde")));

  r = the_seqset->find(tseq("abcde") + dna_A);
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_A, tseq("abcde") + dna_A));

  r = the_seqset->find(tseq("abcde") + dna_C);
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_C));

  r = the_seqset->find(tseq("abcde") + dna_C + tseq("f"));
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_C,
                                   tseq("abcde") + dna_C + tseq("f")));

  r = the_seqset->find(tseq("abcde") + dna_T);
  EXPECT_THAT(get_prefix_read_seqs(the_readmap.get(), r),
              UnorderedElementsAre(tseq("abcde"), tseq("abcde") + dna_T));
}

TEST(readmap, min_max_size) {
  std::vector<std::vector<dna_sequence>> test_reads =  //
      {{tseq("abcde")},
       {tseq("abcde") + dna_A},
       {tseq("abcde") + dna_A},
       {tseq("abcde") + dna_C},
       {tseq("abcde") + dna_C + tseq("f")},
       {tseq("abcde") + dna_T}};
  auto bgfile = biograph_for_reads(test_reads);
  std::shared_ptr<seqset> the_seqset = std::move(bgfile.first);
  std::unique_ptr<readmap> the_readmap = std::move(bgfile.second);

  EXPECT_EQ(the_readmap->min_read_len(), tseq("abcde").size());
  EXPECT_EQ(the_readmap->max_read_len(), (tseq("abcde") + dna_C + tseq("f")).size());
}

TEST(readmap, get_prefix_read_seqs_wild) {
  biograph bg("golden/e_coli_merged.bg");

  std::shared_ptr<readmap> rm1 = bg.open_readmap("e_coli_test");
  std::shared_ptr<readmap> rm2 = bg.open_readmap("test_accession_id");
  std::shared_ptr<seqset> ss = bg.get_seqset();
  seqset_range r = ss->find(dna_sequence("AAAATTACAGAGTACACAACATCCATGAAACGCAT"));
  ASSERT_TRUE(r.valid());

  std::set<int> read_ids;
  for (const auto& read : rm1->get_prefix_reads(r)) {
    read_ids.insert(read.get_read_id());
  }
  EXPECT_THAT(read_ids, IsEmpty());

  read_ids.clear();
  for (const auto& read : rm2->get_prefix_reads(r)) {
    read_ids.insert(read.get_read_id());
  }
  EXPECT_THAT(read_ids, UnorderedElementsAre(123, 124));

  read_ids.clear();
  r = ss->find(dna_sequence("A"));
  for (const auto& read : rm1->get_prefix_reads(r)) {
    read_ids.insert(read.get_read_id());
  }
  EXPECT_THAT(read_ids, IsEmpty());
}

void test_increase_prefix_search_range(seqset_range r, const readmap& rm) {
  SCOPED_TRACE(r.sequence());
  std::set<uint32_t> read_ids;

  for (const auto& read : rm.get_prefix_reads(r)) {
    EXPECT_TRUE(read_ids.insert(read.get_read_id()).second);
  }

  while (r.size() > 0) {
    r = r.truncate(r.size() - 1);
    std::set<uint32_t> new_read_ids;
    for (const auto& read : rm.get_prefix_reads(r)) {
      EXPECT_TRUE(new_read_ids.insert(read.get_read_id()).second);
      // Don't expect any new reads to show up
      EXPECT_TRUE(read_ids.count(read.get_read_id()));
    }

    read_ids = std::move(new_read_ids);
  }
  EXPECT_THAT(read_ids, IsEmpty());
}

// Tests to make sure we don't find any more reads as we increase the search range
TEST(readmap, get_prefix_read_seqs_wild_rand) {
  biograph bg("golden/e_coli_merged.bg");

  std::random_device rand_dev;
  auto seed = rand_dev();
  std::cout << "Rand seed: " << seed << "\n";
  std::mt19937 rand_source(seed);

  std::shared_ptr<seqset> ss = bg.get_seqset();
  std::shared_ptr<readmap> rm;
  if (rand_source()&1) {
    rm = bg.open_readmap("e_coli_test");
  } else {
    rm = bg.open_readmap("e_coli_test");
  }

  uint64_t seqset_id = std::uniform_int_distribution<uint64_t>(0, ss->size() - 1)(rand_source);
  seqset_range r = ss->ctx_entry(seqset_id);

  test_increase_prefix_search_range(r, *rm);
}

// Tests to make sure each read is present in all its suffixes.
TEST(readmap, get_prefix_read_seqs_wild_rand2) {
  biograph bg("golden/e_coli_merged.bg");

  std::random_device rand_dev;
  auto seed = rand_dev();
  std::cout << "Rand seed: " << seed << "\n";
  std::mt19937 rand_source(seed);

  std::shared_ptr<seqset> ss = bg.get_seqset();
  std::shared_ptr<readmap> rm;
  if (rand_source()&1) {
    rm = bg.open_readmap("e_coli_test");
  } else {
    rm = bg.open_readmap("e_coli_test");
  }

  uint32_t read_id = std::uniform_int_distribution<uint32_t>(0, rm->size() - 1)(rand_source);
  int read_len = rm->get_readlength(read_id);
  seqset_range r = ss->ctx_entry(rm->index_to_entry(read_id));

  while (int(r.size()) >= read_len) {
    std::set<uint32_t> read_ids;
    for (const auto& read : rm->get_prefix_reads(r)) {
      EXPECT_TRUE(read_ids.insert(read.get_read_id()).second);
    }
    EXPECT_THAT(read_ids, Contains(read_id));

    test_increase_prefix_search_range(r, *rm);

    r = r.truncate(r.size() - 1);
  }
}

class reads_containing_test : public Test {
 public:
  void SetUp() {
    std::vector<std::vector<dna_sequence>> reads = {
        {tseq("abcde")},        //
        {tseq("abcd")},         //
        {tseq("fazf")},         //
        {tseq("fooaafoo")},     //
        {tseq_rc("fooaafoo")},  //
        {tseq("endswitha")},    //
        {tseq("endswith")},     //
        {tseq("doesnotm*tch")}  //
    };
    auto bgfile = biograph_for_reads(reads);
    m_seqset = std::move(bgfile.first);
    m_readmap = std::move(bgfile.second);
  }

  void search_for(dna_sequence search_seq) {
    CHECK(m_actual.empty());
    seqset_range r = m_seqset->find(search_seq);
    for (std::pair<int /* offset */, readmap::read> offset_and_read :
         m_readmap->get_reads_containing(r)) {
      int offset = offset_and_read.first;
      dna_sequence actual_seq = offset_and_read.second.get_seqset_entry().sequence();
      EXPECT_EQ(actual_seq.subseq(offset, search_seq.size()), search_seq);
      m_actual.emplace_back(std::make_pair(offset, actual_seq));
    }
  }

 protected:
  std::shared_ptr<seqset> m_seqset;
  std::unique_ptr<readmap> m_readmap;

  std::vector<std::pair<int /* offset */, dna_sequence>> m_actual;
};

TEST_F(reads_containing_test, lots) {
  search_for(tseq("a"));

  EXPECT_THAT(m_actual, UnorderedElementsAre(Pair(0, tseq("abcde")), Pair(0, tseq("abcd")),
                                             Pair(tseq("f").size(), tseq("fazf")),
                                             Pair(tseq("foo").size(), tseq("fooaafoo")),
                                             Pair(tseq("foo").size(), tseq("fooaafoo")),
                                             Pair(tseq("fooa").size(), tseq("fooaafoo")),
                                             Pair(tseq("fooa").size(), tseq("fooaafoo")),
                                             Pair(tseq("endswith").size(), tseq("endswitha"))));
}

TEST_F(reads_containing_test, missing) {
  search_for(tseq("A"));

  EXPECT_THAT(m_actual, IsEmpty());
}

TEST_F(reads_containing_test, single_match) {
  search_for(tseq("*"));

  EXPECT_THAT(m_actual, UnorderedElementsAre(Pair(tseq("doesnotm").size(), tseq("doesnotm*tch"))));
}

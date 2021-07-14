#include "base/base.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_mapred/align_seqset.h"
#include "modules/bio_mapred/kmerize_bf.h"
#include "modules/bio_mapred/make_readmap.h"
#include "modules/bio_mapred/mem_seqset.h"
#include "modules/bio_mapred/prefix_flatten.h"
#include "modules/bio_mapred/sort_expand.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/base_chunker.h"
#include "modules/mapred/kv_hold.h"
#include "modules/mapred/resource_manager.h"
#include "modules/mapred/task_mgr.h"
#include "modules/pipeline/primitives.h"
#include "modules/test/build_ref.h"
#include "modules/test/test_utils.h"

#include <bitset>
#include <memory>
#include <random>

#include <gtest/gtest.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

static std::once_flag once_seed_flag;

static void find_inner(std::vector<int>& cov, const std::string& seq, const std::string& read) {
  size_t pos = seq.find(read);
  while (pos != std::string::npos) {
    for (size_t i = 0; i < read.size(); i++) {
      cov[pos + i]++;
    }
    pos = seq.find(read, pos + 1);
  }
}

static void find_cr(std::vector<int>& cov, const std::string& seq, const corrected_reads& crs) {
  for (size_t i = 0; i < crs.size(); i++) {
    std::string fwd = crs[i].corrected.as_string();
    std::string rev = crs[i].corrected.rev_comp().as_string();
    find_inner(cov, seq, fwd);
    find_inner(cov, seq, rev);
  }
}

// NOTE: Not thread safe!
// Given a dna_base, return a different randomly generated base.
static dna_base random_snp(dna_base the_base) {
  static std::random_device entropy_for_seed;
  static unsigned int seed = entropy_for_seed();
  std::call_once(once_seed_flag, []() { SPLOG("Random SNP seed = %u", seed); });

  static std::default_random_engine the_generator(seed);
  static std::uniform_int_distribution<int> random_snp_picker(0, 2);
  auto random_snp = random_snp_picker(the_generator);
  return dna_base((int(the_base) + 1 + random_snp) % 4);
}

static manifest make_cr_manifest(const std::vector<corrected_reads>& the_corrected_reads,
                                 size_t goal_size) {
  manifest corrected_reads_manifest;
  output_stream_params the_output_params;
  the_output_params.goal_size = goal_size;
  auto manifest_kv_sink_ptr = the_output_params.build(make_path("corrected_reads_files"),
                                                      "corrected_reads", corrected_reads_manifest);

  CHECK(manifest_kv_sink_ptr.get());
  int i = 0;
  for (const auto& some_corrected_reads : the_corrected_reads) {
    manifest_kv_sink_ptr->write_msgpack(std::to_string(i++), some_corrected_reads);
  }
  manifest_kv_sink_ptr->close();

  return corrected_reads_manifest;
}

// Function assumes unpaired reads only and constant size reads."
void write_reads_as_fastq(const std::vector<corrected_reads>& the_reads) {
  path fastq_path(make_path("fastq_reads"));
  file_writer fastq_writer(fastq_path.bare_path());

  auto read_id = 0;
  std::string fake_quality;
  for (const auto& the_corrected_reads : the_reads) {
    std::string read_id_string(std::to_string(read_id++));
    fastq_writer.write("@", 1);
    fastq_writer.write(read_id_string.data(), read_id_string.size());
    fastq_writer.write("\n", 1);
    fastq_writer.write(the_corrected_reads[0].corrected.as_string().data(),
                       the_corrected_reads[0].corrected.size());
    fastq_writer.write("\n+", 2);
    fastq_writer.write(read_id_string.data(), read_id_string.size());
    fastq_writer.write("\n", 1);
    fake_quality = std::string(the_reads[0][0].corrected.size(), 'E');
    fastq_writer.write(fake_quality.data(), fake_quality.size());
    fastq_writer.write("\n", 1);
  }

  fastq_writer.close();
}

class seqset_test : public testing::Test {
 public:
  static void SetUpTestCase() { add_primitives(); }
};

TEST_F(seqset_test, mem) {
  perform_build_ref("e_coli", "datasets/fasta/e_coli_k12.ASM584v1.fasta");
  task_mgr_local local_task_manager;

  manifest correct_reads_manifest;
  correct_reads_manifest.add(
      file_info(path("datasets/reads/e_coli_corrected_reads.kvp"), 9731499, 53238), 0);
  correct_reads_manifest.metadata().set(meta::ns::readonly, "corrected_read_bases", 9731499);

  SPLOG("Generating SEQSET");
  auto seqset_task = make_unique<mem_seqset_task>();
  seqset_task->input = correct_reads_manifest;
  seqset_task->num_threads = 32;
  seqset_task->max_mem = 16;
  seqset_task->run_tests = true;
  seqset_task->ref_name = "e_coli";
  manifest seqset_manifest;
  path seqset_out_path(make_path("mem_seqset_task"));
  local_task_manager.run_task(seqset_manifest, seqset_out_path, std::move(seqset_task));
  std::string out_name = seqset_manifest.begin()->file.bare_path();

  auto the_seqset = std::make_shared<seqset>(out_name);

  seqset_range entire_graph{the_seqset.get()};
  ASSERT_TRUE(entire_graph.valid());

  seqset_range first_context = the_seqset->ctx_begin();
  ASSERT_TRUE(first_context.valid());

  dna_sequence probe_read_sequence = std::string(
      "TCAGACTTGATACATTTTAGTTACATATATTTTCTTATTTTATGC"
      "GGAAAATGCTATATGGAAATGTAGTAATTATATACATCTTATCGAAAGTGATTTT");
  seqset_range probe_context = the_seqset->find(probe_read_sequence);
  SPLOG("Probe context size = %d, begin = %lu, end = %lu", probe_context.size(),
        probe_context.begin(), probe_context.end());
  ASSERT_TRUE(probe_context.valid());
  ASSERT_EQ(probe_context.end() - probe_context.begin(), 1);
  ASSERT_EQ(probe_context.size(), probe_read_sequence.size());

  uint64_t previous_begin = 0;
  uint64_t previous_end = the_seqset->size();
  for (unsigned i = 0; i < probe_read_sequence.size(); i++) {
    seqset_range prefix_probe_context = the_seqset->find(probe_read_sequence.subseq(0, i));
    ASSERT_LE(previous_begin, prefix_probe_context.begin());
    ASSERT_GE(previous_end, prefix_probe_context.end());
    ASSERT_TRUE(prefix_probe_context.valid());
    previous_begin = prefix_probe_context.begin();
    previous_end = prefix_probe_context.end();
  }

  uint64_t previous_size = 1;
  for (unsigned i = 0; i < probe_read_sequence.size(); i++) {
    seqset_range suffix_probe_context = the_seqset->find(probe_read_sequence.subseq(i, 100 - i));
    ASSERT_LE(previous_size, suffix_probe_context.end() - suffix_probe_context.begin());
    ASSERT_TRUE(suffix_probe_context.valid());
    previous_size = suffix_probe_context.end() - suffix_probe_context.begin();
  }

  std::vector<std::string> not_in_seqset_vector = {
      {"ATACTAGACAGTAAATAAAATTTTCCTTTGTTCCAGAAGGAGGTACTG"
       "GTTTTCTATTCCAAGGGTGTTTTCTATACAAACATGCTTGAAAATAATCATT"},
      {"GCTTCAGTCTCCCAAGTATTTGGAACTATAAGGTGAACACCACCATACCTGGC"
       "TATTTTTGTTACTTTTTATTTTGTAGAGATGGGGTCTTGCTGTGTTG"},
      {"TCTTGGAGAGGGCCAAGACACTACATGGCCCAGAAGATCACAGTCAGGAGAAA"
       "TACCTGAGCATCTCACAGGACAGATCTGGTGGAAATACCGCTCTGCT"},
      {"GGACGAGCCGCCCCGGCGGTGAACGGGGAGGAGGCGGGAACCGAAGAAGCGG"
       "GGGCGCCGGCCGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCC"},
      {"ACTTATTTAATTTCATTAAAAAACTATCTGAATGCCTCCTTTGTGCAAGATA"
       "TTTGCAAGACAGTGCAAATTGATACAGAAGCTAGTAACACATGGCCCT"}};

  for (auto not_in_seqset : not_in_seqset_vector) {
    seqset_range not_found_context = the_seqset->find(not_in_seqset);
    ASSERT_FALSE(not_found_context.valid());
  }
}

void test_seqset_construct(std::string refname, bool is_paired, unsigned int seed) {
  std::mt19937 the_generator(seed);
  reference hiv_ref(refname);
  dna_sequence modified_hiv_ref_seq(hiv_ref.get_dna(0), hiv_ref.get_dna(0) + hiv_ref.size());

  constexpr int coverage = 50;
  constexpr unsigned read_size = 100;
  constexpr int read_size_delta = 20;  // Read size will vary up or down by this amount.
  constexpr int snp_count = 100;
  constexpr size_t goal_size = 1000;
  // constexpr size_t goal_size = 64 * 1024 * 1024; // Set to 64MB to make the test single threaded
  constexpr int kmer_size = 8;
  std::uniform_int_distribution<int> snp_location_picker(read_size, hiv_ref.size() - read_size - 1);
  std::uniform_int_distribution<int> direction_picker(0, 1);  // 0 = forward, 1 = reverse
  std::uniform_int_distribution<int> delta_picker(-read_size_delta, read_size_delta);

  SPLOG("Testing with pairing %s", is_paired ? "on" : "off");

  for (auto i = 0; i < snp_count; i++) {
    auto snp_location = snp_location_picker(the_generator);
    dna_base snp_base = random_snp(modified_hiv_ref_seq[snp_location]);
    SPLOG("Changing base %d from %c to %c", snp_location, char(modified_hiv_ref_seq[snp_location]),
          char(snp_base));
    modified_hiv_ref_seq[snp_location] = snp_base;
  }

  auto read_count = coverage * modified_hiv_ref_seq.size() / read_size;
  SPLOG(
      "Ref size = %lu, Coverage = %d, Read size = %u, Read size delta = %d, SNP count = %d, Read "
      "count = %lu",
      modified_hiv_ref_seq.size(), coverage, read_size, read_size_delta, snp_count, read_count);
  std::vector<corrected_reads> the_corrected_reads;
  the_corrected_reads.reserve(read_count);

  auto total_bases = 0;
  auto ref_iter = modified_hiv_ref_seq.begin();
  for (auto i = 0UL; i < read_count; i++) {
    int read_location = snp_location_picker(the_generator);
    int read_direction = direction_picker(the_generator);
    int actual_read_size = read_size + delta_picker(the_generator);
    if (is_paired && (i % 2) == 0) {
      the_corrected_reads.emplace_back(std::vector<corrected_read>(2));
    } else if (is_paired && (i % 2) == 1) {
      // Do nothing, use second element of previously emplaced pair.
    } else if (!is_paired) {
      the_corrected_reads.emplace_back(std::vector<corrected_read>(1));
    } else {
      throw io_exception("Programmer error: this code should be unreachable.");
    }
    if (read_direction == 0) {
      if (read_location + actual_read_size > int(modified_hiv_ref_seq.size())) {
        SPLOG(
            "Skipping bad SNP attempted on iteration %lu at ref location %d, actual read length = "
            "%u",
            i, read_location, actual_read_size);
        the_corrected_reads.pop_back();
        i--;
        continue;
      }
      the_corrected_reads.back()[is_paired ? i % 2 : 0].corrected =
          dna_sequence(ref_iter + read_location, ref_iter + read_location + actual_read_size);
    }

    if (read_direction == 1) {
      if (read_location + 1 < actual_read_size) {
        SPLOG(
            "Skipping reversed bad SNP attempted on iteration %lu at ref location %d, actual read "
            "length = %u",
            i, read_location, actual_read_size);
        the_corrected_reads.pop_back();
        i--;
        continue;
      }
      the_corrected_reads.back()[is_paired ? i % 2 : 0].corrected =
          dna_sequence(ref_iter + read_location + 1 - actual_read_size,
                       ref_iter + read_location + 1)
              .rev_comp();
    }
    ASSERT_GT(the_corrected_reads.back()[is_paired ? i % 2 : 0].corrected.size(), 0);
    total_bases += actual_read_size;
  }
  manifest corrected_reads_manifest = make_cr_manifest(the_corrected_reads, goal_size);
  corrected_reads_manifest.metadata().set(meta::ns::readonly, "corrected_read_bases", total_bases);
  SPLOG("Corrected reads manifest has %lu file_infos", corrected_reads_manifest.count_file_infos());

  write_reads_as_fastq(the_corrected_reads);

  auto build_seqset_task = make_unique<mem_seqset_task>();
  build_seqset_task->input = corrected_reads_manifest;
  build_seqset_task->is_paired = is_paired;
  build_seqset_task->ref_name = "hiv";
  build_seqset_task->num_threads = 32;
  build_seqset_task->max_mem = 8;
  build_seqset_task->run_tests = true;

  task_mgr_local local_task_manager;
  manifest seqset_manifest;
  path seqset_out_path(make_path("construct_seqset_task"));
  local_task_manager.run_task(seqset_manifest, seqset_out_path, std::move(build_seqset_task));

  auto the_seqset = std::make_shared<seqset>(seqset_manifest.begin()->file.bare_path());
  seqset_range entire_graph{the_seqset.get()};
  ASSERT_TRUE(entire_graph.valid());

  SPLOG("Making the readmap")
  std::string readmap_filename =
      is_paired ? make_path("paired_readmap") : make_path("unpaired_readmap");
  make_readmap::do_make(readmap_filename, *the_seqset, corrected_reads_manifest, is_paired,
                        std::numeric_limits<uint8_t>::max());

  SPLOG("Loading readmap");
  readmap the_readmap(the_seqset, readmap_filename);

  SPLOG("Looking for existing reads.");
  for (const auto& corrected_reads : the_corrected_reads) {
    seqset_range read_range = the_seqset->find(corrected_reads[0].corrected);
    ASSERT_TRUE(read_range.valid());

    auto entry_id = read_range.begin();  // Used in readmap tests below
    dna_slice read_slice = dna_slice(corrected_reads[0].corrected);

    dna_slice test_slice = read_slice;
    while (test_slice.size() > 10) {  // Prefixes
      seqset_range old_range = read_range;
      read_range = the_seqset->find(test_slice);
      ASSERT_TRUE(read_range.valid());
      ASSERT_LE(read_range.begin(), old_range.begin());
      ASSERT_GE(read_range.end(), old_range.end());
      ASSERT_TRUE(the_seqset->find(test_slice.rev_comp()).valid());
      test_slice = dna_slice(test_slice.begin(), test_slice.size() - 1);
    }

    test_slice = read_slice;
    while (test_slice.size() > 10) {  // Suffixes
      read_range = the_seqset->find(test_slice);
      ASSERT_TRUE(read_range.valid());
      ASSERT_TRUE(the_seqset->find(test_slice.rev_comp()).valid());
      test_slice = dna_slice(test_slice.begin() + 1, test_slice.size() - 1);
    }

    test_slice = read_slice;
    while (test_slice.size() > 10) {  // Midfixes
      read_range = the_seqset->find(test_slice);
      ASSERT_TRUE(read_range.valid());
      ASSERT_TRUE(the_seqset->find(test_slice.rev_comp()).valid());
      test_slice = dna_slice(test_slice.begin() + 1, test_slice.size() - 2);
    }
    // Test readmap matepair information and is_forward
    //
    auto indexes = the_readmap.entry_to_index(entry_id);
    if (indexes.first == indexes.second) {  // single entry to index
      ASSERT_EQ(the_readmap.has_mate(indexes.first), is_paired);
      ASSERT_TRUE(the_readmap.get_is_forward(indexes.first));

      // Make sure the read's reverse complement isn't is_forward...
      // That is if it returns a single entry (there's no forward revcmp seq)
      auto rc_entry_id = the_seqset->find(corrected_reads[0].corrected.rev_comp());
      auto rc_index = the_readmap.entry_to_index(rc_entry_id.begin());
      if (rc_index.first == rc_index.second) {
        ASSERT_FALSE(the_readmap.get_is_forward(rc_index.first));
      }

      // Check the mate
      if (is_paired) {
        auto mate_entry = the_seqset->find(corrected_reads[1].corrected).begin();
        auto mate_indexes = the_readmap.entry_to_index(mate_entry);
        if (mate_indexes.first == mate_indexes.second) {  // single entry to index
          // r1 points to r2
          ASSERT_EQ(the_readmap.get_mate(indexes.first), mate_indexes.first);
          // r2 points to r2
          ASSERT_EQ(the_readmap.get_mate(mate_indexes.first), indexes.first);
          // r1 can recreate correct entry_id
          ASSERT_EQ(the_readmap.index_to_entry(indexes.first), mate_entry);
          // r2 can recreate correct entry_id
          ASSERT_EQ(the_readmap.index_to_entry(mate_indexes.first), entry_id);
        }
      }
    }  // else { Need to figure out which one I'm looking at }

    // Test 'find_near' code
    // First copy the original
    read_range = the_seqset->find(corrected_reads[0].corrected);
    dna_sequence mread = corrected_reads[0].corrected;
    // Now add two random SNPS
    mread[random() % mread.size()] = dna_base(int(random() % 4));
    mread[random() % mread.size()] = dna_base(int(random() % 4));
    std::vector<seqset_range> out;
    bool r = the_seqset->find_near(out, mread, 2, 1000);
    ASSERT_TRUE(r);
    bool got_it = false;
    for (size_t i = 0; i < out.size(); i++) {
      if (out[i] == read_range) {
        got_it = true;
      }
    }
    ASSERT_TRUE(got_it);
  }

  SPLOG("Testing coverage");
  for (size_t i = 0; i < 100; i++) {
    size_t start = random() % (modified_hiv_ref_seq.size() - 700);
    size_t end = start + 200 + random() % 500;
    dna_slice x(modified_hiv_ref_seq.begin() + start, modified_hiv_ref_seq.begin() + end);
    std::vector<int> cov = the_readmap.approx_coverage(x);
    std::vector<int> rcov(x.size());
    std::string xstr = x.as_string();
    for (size_t j = 0; j < the_corrected_reads.size(); j++) {
      find_cr(rcov, xstr, the_corrected_reads[j]);
    }
    bool err = false;
    for (size_t j = 0; j < cov.size(); j++) {
      if (cov[i] != rcov[i]) {
        err = true;
      }
    }
    if (!err) continue;
    SPLOG("%s", x.as_string().c_str());
    std::string depths = "[";
    for (size_t j = 0; j < cov.size(); j++) {
      depths += printstring("%d, ", cov[j]);
    }
    depths = depths.substr(0, depths.size() - 2) + "]";
    SPLOG("%s", x.as_string().c_str());
    std::string depths2 = "[";
    for (size_t j = 0; j < cov.size(); j++) {
      depths2 += printstring("%d, ", rcov[j]);
    }
    depths2 = depths2.substr(0, depths2.size() - 2) + "]";
    SPLOG("%s", depths.c_str());
    SPLOG("%s", depths2.c_str());
    // This one fails
    ASSERT_TRUE(!err);
  }

  std::bitset<kmer_t(1) << (kmer_size * 2)> kmer_bitset;
  for (auto ref_seq_iter = modified_hiv_ref_seq.begin();
       ref_seq_iter < modified_hiv_ref_seq.end() - kmer_size + 1; ref_seq_iter++) {
    kmer_bitset.set(dna_slice(ref_seq_iter, kmer_size).as_kmer());
    kmer_bitset.set(dna_slice(ref_seq_iter, kmer_size).rev_comp().as_kmer());
  }
  SPLOG("%lu kmers of %lu found in reference", kmer_bitset.count(), kmer_bitset.size());

  for (auto i = 0UL; i < kmer_bitset.size(); i++) {
    if (!kmer_bitset.test(i)) {
      auto kmer_range = the_seqset->find(dna_sequence(kmer_t(i), kmer_size));
      if (kmer_range.valid()) {
        SPLOG("i = %lu, Kmer = 0x%lX, DNA = %s", i, i,
              dna_sequence(kmer_t(i), kmer_size).as_string().c_str());
        SPLOG("Begin = %lu, end = %lu, size = %d, sequence = %s", kmer_range.begin(),
              kmer_range.end(), kmer_range.size(), kmer_range.sequence().as_string().c_str());
        // This one fails
        ASSERT_FALSE(true);
      }
    }
  }
}

TEST_F(seqset_test, construct_unpaired) {
  perform_build_ref("hiv", "datasets/hiv/ref/hiv-1-NC_001802.1.fa");

  std::random_device seed_entropy;
  unsigned int seed = seed_entropy();
  // unsigned int seed = 1955666968;
  SPLOG("SNP location seed = %u", seed);

  test_seqset_construct("hiv", false, seed);
}

TEST_F(seqset_test, construct_paired) {
  perform_build_ref("hiv", "datasets/hiv/ref/hiv-1-NC_001802.1.fa");

  std::random_device seed_entropy;
  unsigned int seed = seed_entropy();
  // unsigned int seed = 1955666968;
  SPLOG("SNP location seed = %u", seed);

  test_seqset_construct("hiv", true, seed);
}

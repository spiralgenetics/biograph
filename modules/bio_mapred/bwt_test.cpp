#include <memory>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <gtest/gtest.h>

#include "modules/bio_base/bwt_file.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_mapred/kmerize_bf.h"
#include "modules/bio_mapred/make_bwt.h"
#include "modules/bio_mapred/prefix_flatten.h"
#include "modules/bio_mapred/sort_expand.h"
#include "modules/pipeline/primitives.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/base_chunker.h"
#include "modules/mapred/kv_hold.h"
#include "modules/mapred/resource_manager.h"
#include "modules/mapred/task_mgr.h"

#include "modules/test/build_ref.h"
#include "modules/test/test_utils.h"

// simple.fasta looks like this:
//
//           1         2
// 0123456789012345678901234567
// ATTGCTC
//        AGAGCTCACTG
//                   TGACTCACTG
// SC:1   SC:2       SC:N

void run_queries(const std::string& bwt) {
	bwt_file my_bwt(bwt);

	size_t matches;

	// exact match
	matches = my_bwt.bwt().find(dna_sequence("GACTCACT")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found GACTCACT: %lu", matches);

	// multiple match
	matches = my_bwt.bwt().find(dna_sequence("A")).matches();
	ASSERT_EQ(matches, 6);
	SPLOG("Found A: %lu", matches);
	matches = my_bwt.bwt().find(dna_sequence("C")).matches();
	ASSERT_EQ(matches, 8);
	SPLOG("Found C: %lu", matches);
	matches = my_bwt.bwt().find(dna_sequence("G")).matches();
	ASSERT_EQ(matches, 6);
	SPLOG("Found G: %lu", matches);
	matches = my_bwt.bwt().find(dna_sequence("T")).matches();
	ASSERT_EQ(matches, 8);
	SPLOG("Found T: %lu", matches);

	// Where are all of the Ts?
	std::vector<size_t>t_positions = {1, 2, 5, 12, 16, 18, 22, 26};
	auto t_bwt = my_bwt.bwt().find(dna_sequence("T"));
	for(size_t i = 0; i < matches; i++) {
		auto it = find(t_positions.begin(), t_positions.end(), t_bwt.get_match(i));
		ASSERT_NE(it, t_positions.end());
		SPLOG("T found at: %lu", *it);
	  	t_positions.erase(it);
	}
	// Did we find them all?
	ASSERT_EQ(t_positions.size(), 0);

	// mismatch
	matches = my_bwt.bwt().find(dna_sequence("GGGGG")).matches();
	ASSERT_EQ(matches, 0);
	SPLOG("Found GGGGG: %lu", matches);

	// match beginning of each sc
	matches = my_bwt.bwt().find(dna_sequence("ATT")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found ATT: %lu", matches);
	matches = my_bwt.bwt().find(dna_sequence("AGA")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found AGA: %lu", matches);
	matches = my_bwt.bwt().find(dna_sequence("TGA")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found TGA: %lu", matches);

	// match end of each sc
	matches = my_bwt.bwt().find(dna_sequence("TGCTC")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found TGCTC: %lu", matches);
	matches = my_bwt.bwt().find(dna_sequence("GCTCACTG")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found GCTCACTG: %lu", matches);
	matches = my_bwt.bwt().find(dna_sequence("ACTCACTG")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found ACTCACTG: %lu", matches);

	// match whole sc
	matches = my_bwt.bwt().find(dna_sequence("TGACTCACTG")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found TGACTCACTG: %lu", matches);

	// matches shouldn't cross supercontig boundaries
	matches = my_bwt.bwt().find(dna_sequence("TCA")).matches();
	ASSERT_EQ(matches, 2);
	SPLOG("Found TCA: %lu", matches);
}

void run_queries_reuse_context(std::string the_bwt) {
	bwt_file my_bwt(the_bwt);

	size_t matches;
	bwt_range bwt = my_bwt.bwt();

	// exact match
	matches = bwt.find(dna_sequence("GACTCACT")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found GACTCACT: %lu", matches);

	// multiple match
	matches = bwt.find(dna_sequence("A")).matches();
	ASSERT_EQ(matches, 6);
	SPLOG("Found A: %lu", matches);
	matches = bwt.find(dna_sequence("C")).matches();
	ASSERT_EQ(matches, 8);
	SPLOG("Found C: %lu", matches);
	matches = bwt.find(dna_sequence("G")).matches();
	ASSERT_EQ(matches, 6);
	SPLOG("Found G: %lu", matches);
	matches = bwt.find(dna_sequence("T")).matches();
	ASSERT_EQ(matches, 8);
	SPLOG("Found T: %lu", matches);

	// Where are all of the Ts?
	std::vector<size_t>t_positions = {1, 2, 5, 12, 16, 18, 22, 26};
	auto t_bwt = bwt.find(dna_sequence("T"));
	for(size_t i = 0; i < matches; i++) {
		auto it = find(t_positions.begin(), t_positions.end(), t_bwt.get_match(i));
		ASSERT_NE(it, t_positions.end());
		SPLOG("T found at: %lu", *it);
	  	t_positions.erase(it);
	}
	// Did we find them all?
	ASSERT_EQ(t_positions.size(), 0);

	// mismatch
	matches = bwt.find(dna_sequence("GGGGG")).matches();
	ASSERT_EQ(matches, 0);
	SPLOG("Found GGGGG: %lu", matches);

	// match beginning of each sc
	matches = bwt.find(dna_sequence("ATT")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found ATT: %lu", matches);
	matches = bwt.find(dna_sequence("AGA")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found AGA: %lu", matches);
	matches = bwt.find(dna_sequence("TGA")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found TGA: %lu", matches);

	// match end of each sc
	matches = bwt.find(dna_sequence("TGCTC")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found TGCTC: %lu", matches);
	matches = bwt.find(dna_sequence("GCTCACTG")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found GCTCACTG: %lu", matches);
	matches = bwt.find(dna_sequence("ACTCACTG")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found ACTCACTG: %lu", matches);

	// match whole sc
	matches = bwt.find(dna_sequence("TGACTCACTG")).matches();
	ASSERT_EQ(matches, 1);
	SPLOG("Found TGACTCACTG: %lu", matches);

	// matches shouldn't cross supercontig boundaries
	matches = bwt.find(dna_sequence("TCA")).matches();
	ASSERT_EQ(matches, 2);
	SPLOG("Found TCA: %lu", matches);
}

TEST(bwt_test, simple)
{
	perform_build_ref("simple", "golden/ftest/bwt/simple.fasta");

	// built automatically by reference build
	run_queries(CONF_S(reference_path) + "/simple/reference.bwt");
	run_queries_reuse_context(CONF_S(reference_path) + "/simple/reference.bwt");
}

TEST(bwt_test, builtin)
{
	perform_build_ref("simple", "golden/ftest/bwt/simple.fasta");
	task_mgr_local local_task_manager;

	// use the bwt built into the reference
	reference ref("simple");

	bwt_range bwt = ref.get_bwt();
	size_t matches;
	matches = bwt.find(dna_sequence("A")).matches();
	ASSERT_EQ(matches, 6);
	SPLOG("Found A: %lu", matches);
	matches = bwt.find(dna_sequence("C")).matches();
	ASSERT_EQ(matches, 8);
	SPLOG("Found C: %lu", matches);
	matches = bwt.find(dna_sequence("G")).matches();
	ASSERT_EQ(matches, 6);
	SPLOG("Found G: %lu", matches);
	matches = bwt.find(dna_sequence("T")).matches();
	ASSERT_EQ(matches, 8);
	SPLOG("Found T: %lu", matches);
}

TEST(bwt_test, small_century)
{
	perform_build_ref("simple", "golden/ftest/bwt/simple.fasta");
	task_mgr_local local_task_manager;

	SPLOG("Generating BWT");
	std::string out = CONF_S(reference_path) + "/simple.bwt";
	auto bwt_task = make_unique<make_bwt_task>();
	bwt_task->input_ref = CONF_S(reference_path) + "/simple/reference.ref";
	bwt_task->output_bwt = out;
	bwt_task->cent_mod = 3;

	path bwt_out_path(make_path("make_bwt_task"));
	local_task_manager.run_task(out, bwt_out_path, std::move(bwt_task));

	run_queries(out);
	run_queries_reuse_context(out);

}

TEST(bwt_test, kmer)
{
	// Negative + positive kmer test.

	constexpr int kmer_size = 8;

	// Note that HIV only has one supercontig. If you want to use a different
	// organism with multiple chromosomes, you'll need to respect extent boundaries
	// when kmerizing.

	perform_build_ref("hiv", "datasets/hiv/ref/hiv-1-NC_001802.1.fa");
	reference ref("hiv");

	bwt_range bwt = ref.get_bwt();
	dna_sequence ref_seq(ref.get_dna(0), ref.get_dna(0) + ref.size());

	std::bitset<kmer_t(1) << (kmer_size * 2)> kmer_bitset;
	for (auto ref_seq_iter = ref_seq.begin()
		; ref_seq_iter < ref_seq.end() - kmer_size + 1
		; ref_seq_iter++) {
		kmer_bitset.set(dna_slice(ref_seq_iter, kmer_size).as_kmer());
	}
	SPLOG("%lu kmers of %lu found in reference", kmer_bitset.count(), kmer_bitset.size());

	for (auto i = 0UL; i < kmer_bitset.size(); i++) {
		// negative test
		if (! kmer_bitset.test(i)) {
			auto kmer_range = bwt.find(dna_sequence(kmer_t(i), kmer_size));
			if (kmer_range.valid()) {
				SPLOG("Found invalid kmer: i = %lu, Kmer = 0x%lX, DNA = %s", i, i, dna_sequence(kmer_t(i), kmer_size).as_string().c_str());
				SPLOG("Begin = %lu, end = %lu, matches = %lu, sequence = %s"
					, kmer_range.begin(), kmer_range.end(), kmer_range.matches(), dna_sequence(kmer_t(i), kmer_size).as_string().c_str());
				ASSERT_FALSE(true);
			}
		// positive test
		} else {
			auto kmer_range = bwt.find(dna_sequence(kmer_t(i), kmer_size));
			if (! kmer_range.valid()) {
				SPLOG("Missed valid kmer: i = %lu, Kmer = 0x%lX, DNA = %s", i, i, dna_sequence(kmer_t(i), kmer_size).as_string().c_str());
				SPLOG("Begin = %lu, end = %lu, matches = %lu, sequence = %s"
					, kmer_range.begin(), kmer_range.end(), kmer_range.matches(), dna_sequence(kmer_t(i), kmer_size).as_string().c_str());
				ASSERT_FALSE(true);
			}
		}
	}
}

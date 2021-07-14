#pragma once

#include "modules/bio_base/corrected_read.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/io/transfer_object.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/mapred/sorter.h"
#include "modules/io/log.h"

// NOTE: This is for the correct_reads_only step. For correct_reads, see anchored_assembly.h
struct correct_reads_params
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(kmer_db);
		FIELD(max_quality_cost, TF_STRICT);
		FIELD(min_base_quality, TF_STRICT);
		FIELD(trim);
		FIELD(trace);
		FIELD(skip_snps);
		FIELD(exact);
        FIELD(trim_after_portion);
		FIELD(frc_max_corrections);
		FIELD(frc_min_good_run);
	}

	double max_quality_cost;
	double min_base_quality;
	bool trace = false;
	std::string kmer_db;
	size_t trim;
	bool skip_snps = false;
	bool exact = false;

	// When used with exact, trims any reads that fail kmer
	// verification after this portion of the read instead of
	// discarding them entirely.
	double trim_after_portion = 1;

  unsigned frc_max_corrections = 2;
  unsigned frc_min_good_run = 2;

	void validate() {
		SPLOG_P(LOG_DEBUG, "correct_reads_params::validate> min_base_quality: %0.0f, max_quality_cost: %0.0f, trim: %lu, trace: %s, skip_snps: %s, exact: %s",
			min_base_quality,
			max_quality_cost,
			trim,
			trace ? "true" : "false",
			skip_snps ? "true" : "false",
			exact ? "true" : "false"
		);
	}
};

struct read_correction_stats
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(corrected_read_count, 0ULL);
		FIELD(corrected_read_bases, 0ULL);
		FIELD(corrected_base_dist);
		FIELD(failed_correction_count, 0ULL);
	}

	uint64_t corrected_read_count = 0;
	uint64_t corrected_read_bases = 0;
	// Holds the number of reads corrected with the index number of bases.
	// E.g. element 3 contains the number of reads with three bases corrected.
	std::vector<uint64_t> corrected_base_dist;
	uint64_t failed_correction_count = 0;
};

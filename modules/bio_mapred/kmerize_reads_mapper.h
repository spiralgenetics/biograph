#pragma once

#include "modules/io/transfer_object.h"
#include "modules/mapred/mapper.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"

#include <array>

struct kmerize_reads_params
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(kmer_size, TF_STRICT);
		FIELD(overrep, size_t(1024*1024*1024));
		FIELD(prior_count, size_t(5));
		FIELD(skew_cutoff, float(0.0));
		FIELD(use_score, false);
		FIELD(trim);
	}

	size_t kmer_size;
        size_t overrep = 1024*1024*1024;
        size_t prior_count = 5;
        float skew_cutoff = 0.0;
	bool use_score;
	size_t trim;

	void validate();
};

class kmerize_reads_mapper : public typed_mapper<kmerize_reads_mapper,
	read_id, unaligned_reads,
	kmer_t, kcount_pair>
{
public:
	kmerize_reads_mapper(const std::string& params);

	void typed_map(const read_id& key, const unaligned_reads& value);
	task_requirements get_requirements() override;

	static std::array<int64_t, 127> mg_log_lookup_table;

private:
	void map_one_read(const read_id& read_id, const unaligned_read& r);
	// Calculate the quality of a kmer with the given start position in a sequence.
	int64_t check_qual(const read_id& read_id, const std::string& qual, size_t start);
	// Given a kmer, calculate its quality by dropping the first base and appending a new base.
	int64_t rotate_qual(const read_id& read_id, const std::string& qual, size_t start, int64_t log_prob);

private:
	kmerize_reads_params m_params;
};

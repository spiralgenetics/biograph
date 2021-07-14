#pragma once

#include "modules/io/transfer_object.h"
#include "modules/mapred/mapper.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"

class kmer_hist_mapper : public typed_mapper<kmer_hist_mapper,
	kmer_t, kcount_pair,
	uint64_t, uint64_t>
{
public:
	inline kmer_hist_mapper(const std::string& params) {}

	inline void typed_map(kmer_t key, const kcount_pair& value) {
		output(value.fwd + value.rev, 1);
	}
};


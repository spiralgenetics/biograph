#pragma once

#include "modules/bio_mapred/kmer_set.h"
#include "modules/io/transfer_object.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/mapred/sorter.h"

// Returns the number of bases of "read" that match kmers in
// the kmer set.
unsigned verify_kmers(const dna_sequence& read, const kmer_set& kmers);

double align_kmer(std::vector<kmer_t>& path, const dna_sequence& read, const std::string& qual, 
			const kmer_set& kmers, double min_base_quality, double max_cost);

dna_sequence get_corrected(const std::vector<kmer_t>& path, size_t kmer_size);


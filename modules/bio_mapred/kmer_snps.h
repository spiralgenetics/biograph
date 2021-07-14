#pragma once

#include "modules/bio_mapred/kmer_set.h"
#include <map>

void kmer_find_snps(const kmer_set& kmers, std::map<kmer_t, kmer_t>& out, size_t max_memory, size_t num_threads);


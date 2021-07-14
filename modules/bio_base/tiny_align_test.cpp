
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"
#include <gtest/gtest.h>
#include <bitset>
#include <unordered_set>
#include <deque>

constexpr kmer_t fives = 0x5555555555555555l;

#include <boost/multiprecision/cpp_int.hpp>

typedef boost::multiprecision::uint128_t big_kmer_t;
typedef std::vector<big_kmer_t> big_kmers_t;

/*
class patterns
{
public:
	patterns(const std::vector<dna_sequnce> patterns)
	{
		m_fives.resize(patterns.size());
		m_matches.resize(patterns.size());
		m_masks.resize(patterns.size());
		for(size_t i = 0; i < patterns.size(); i++)
		{
			for(size_t 

	big_kmers_t apply_base(const big_kmers_t& errs, int base);
	big_kmers_t m_fives;
	big_kmers_t m_matches;
	big_kmers_t m_masks;
};

big_kmers_t apply_base(const big_kmers_t& errs, const big_kmers_t& patterns, const big_kmers_t& patterns,
*/

kmer_t apply_base(kmer_t errs, kmer_t match, kmer_t mask, kmer_t base) 
{
	errs >>= 2;
	kmer_t base_rep = fives * base;
	kmer_t diff = base_rep ^ match;
	kmer_t diff_bit = (diff | (diff >> 1)) & fives;
	kmer_t saturate = errs & (errs >> 1);
	kmer_t to_add = diff_bit & ~saturate;
	errs += to_add;
	errs &= mask;
	return errs;
}

TEST(tiny_align, test_it) 
{
	dna_sequence match_seq("CTGTCTCTTATACACATCT");
	dna_sequence seek_seq("ACCGTCTGTCTCTTATTACTGTCTCTTATACACATCTGGGTAGA");
	size_t kmer_size = match_seq.size();
	kmer_t mask = (uint64_t(1) << (2*kmer_size)) - 1;
	kmer_t errs = mask;
	kmer_t match = match_seq.as_kmer();
	for(size_t i = 0; i < seek_seq.size(); i++) {
		kmer_t base = (int) seek_seq[i];
		errs = apply_base(errs, match, mask, base);
		printf("%c: %s\n", char(seek_seq[i]),
			std::bitset<64>(errs).to_string().c_str());
	}
}

TEST(tiny_align, test_combin)
{
	dna_sequence match_seq("CTGTCTCTTATACACATCT");
	size_t kmer_size = match_seq.size();
	kmer_t mask = (uint64_t(1) << (2*kmer_size)) - 1;
	kmer_t match = match_seq.as_kmer();
	std::set<kmer_t> found;
	std::deque<kmer_t> to_do;
	found.emplace(mask);
	to_do.push_back(mask);
	while(!to_do.empty()) {
		if (found.size() % 1000 == 0) {
			printf("Found = %d: to_do = %d", int(found.size()), int(to_do.size()));
		}
		kmer_t errs = to_do.front();
		to_do.pop_front();
		for(kmer_t b = 0; b < 4; b++) {
			kmer_t new_err = apply_base(errs, match, mask, b);
			if (found.count(new_err)) {
				continue;
			}
			found.emplace(new_err);
			to_do.push_back(new_err);
		}
	}
	printf("Total found size: %d\n", int(found.size()));	
	for(kmer_t k : found) {
		printf("%s\n", std::bitset<64>(k).to_string().c_str());
	}
}


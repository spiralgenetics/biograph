#include "modules/test/test_utils.h"
#include "modules/bio_base/kmer_counter.h"
#include "modules/io/hash.h"

#include "gtest/gtest.h"

#include <random>

TEST(kmer_counter, basic)
{
	std::mt19937 generator;
	std::uniform_int_distribution<kmer_t> uniform;
	std::vector<kmer_t> samples;
	for (size_t i = 0; i < 50000; i++) {
		samples.push_back(uniform(generator));
	}

	std::map<kmer_t, size_t> counts;
	kmer_counter<basic_hasher<prime_hasher>> kc(100000);
	for (auto kmer : samples) {
		kc.add(kmer, true);
		counts[kmer]++;
	}

	size_t count = 0;
	for (const auto& item : kc) {
		// SPLOG("Key: %ld Value: %u", item.key, item.count);
		count++;
		ASSERT_EQ(counts[item.key], item.fwd_count);
	}

	ASSERT_EQ(counts.size(), count);

	for (auto kmer : samples) {
		kc.add(kmer, true);
		counts[kmer]++;
	}

	count = 0;
	for (const auto& item : kc) {
		count++;
		ASSERT_EQ(counts[item.key], item.fwd_count);
	}

	ASSERT_EQ(counts.size(), count);
}

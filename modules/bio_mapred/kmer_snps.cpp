
#include "modules/bio_mapred/kmer_snps.h"
#include "modules/bio_base/kmer.h"
#include "modules/io/log.h"

typedef std::pair<kmer_t, kmer_t> kmer_match_t;
typedef std::vector<kmer_match_t> kmer_matches_t;

static void exhaustive_search(kmer_matches_t& out, kmer_t k, const kmer_set& ks, uint32_t kmer_size) 
{
	// Modify reverse complement and original simultanously, since bit flip + < is cheaper than RC
	kmer_t kc = rev_comp(k, kmer_size);
	//SPLOG("K  = %s", dna_sequence(k, kmer_size).as_string().c_str());
	//SPLOG("KC = %s", dna_sequence(kc, kmer_size).as_string().c_str());
	// Go over each position
	for(uint32_t i = 0; i < (kc == k ? (kmer_size + 1)/2 : kmer_size); i++) {
		// Go over each modification
		for(uint64_t x = 1; x <= 3; x++) {
			// Modify original + complement
			kmer_t kd = k ^ (x << (2*i));
			kmer_t kcd = kc ^ (x << (2*(kmer_size - i - 1)));
			//SPLOG("   KD  = %s", dna_sequence(kd, kmer_size).as_string().c_str());
			//SPLOG("   KCD = %s", dna_sequence(kcd, kmer_size).as_string().c_str());
			// Pick the canonical one
			kmer_t best = std::min(kd, kcd);
			//SPLOG("   B   = %s", dna_sequence(best, kmer_size).as_string().c_str());
			// Look it up and output if needed
			if (ks.count(best)) {
				//SPLOG("      w00t");
				out.emplace_back(k, best);
				out.emplace_back(best, k);
			}
		}
		
	}
}

void kmer_find_snps(const kmer_set& kmers, std::map<kmer_t, kmer_t>& out, size_t max_memory, size_t num_threads)
{
	SPLOG("Hello, total size = %d", int(kmers.size()));
	kmer_matches_t r;
	uint32_t kmer_size = kmers.kmer_size();
	int count = 0;
	time_t now = time(0);
	for(kmer_t k : kmers) {
		if (time(0) - now > 1) {
			SPLOG("Count = %d", int(count));
			now = time(0);
		}
		exhaustive_search(r, k, kmers, kmer_size);
		count++;
	}
	SPLOG("Done: %d entries", (int) r.size());
	std::sort(r.begin(), r.end());
	for(size_t i = 0; i < r.size(); i++) {
		if (i % 10000 == 0) {
		SPLOG("%s->%s", 
			dna_sequence(r[i].first, kmer_size).as_string().c_str(),
			dna_sequence(r[i].second, kmer_size).as_string().c_str());
		}
	}
	for(size_t i = 0; i < r.size()/2; i++) {
		if (r[2*i] != r[2*i+1]) {
			throw io_exception(printstring("Mismatch: %s->%s vs %s->%s",
				dna_sequence(r[2*i].first, kmer_size).as_string().c_str(),
				dna_sequence(r[2*i].second, kmer_size).as_string().c_str(),
				dna_sequence(r[2*i+1].first, kmer_size).as_string().c_str(),
				dna_sequence(r[2*i+1].second, kmer_size).as_string().c_str()));
		}
	}
}

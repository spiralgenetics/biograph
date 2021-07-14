
#pragma once

#include "modules/bio_base/corrected_read.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/io/transfer_object.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/mapred/mapper.h"

struct filter_reads_params 
{
	TRANSFER_OBJECT 
	{ 
		VERSION(0); 
		FIELD(kmer_db);   // For kmer alignment
		FIELD(hard_filter);  // Delete or mark trace
	}

	std::string kmer_db; 
	bool hard_filter = false;

	void validate() {}
}; 

class filter_reads_mapper: public typed_mapper<filter_reads_mapper,
	std::string, corrected_reads, 
	std::string, corrected_reads>
{
public:
	filter_reads_mapper(const std::string& params);

	void set_watchdog(const std::function<void()>& watchdog) override
	{ 
		m_watchdog = watchdog; 
	}

	void setup() override;
	void typed_map(const std::string& key, const corrected_reads& value);	
	void install_metadata(meta::data& metadata) override;

private:
	bool map_one_read(const std::string& key, corrected_read& r);
	std::function<void()> m_watchdog;
	filter_reads_params m_params;
	std::unique_ptr<const kmer_set> m_kdb;
	// k-mer count in a read not matching DB distribution.  Element 0 is number with
	// 0 mismatches, element 1 is single-mismatch count etc.
	std::vector<unsigned> m_mismatch_count;
	size_t m_tot_mapped;
	size_t m_tot_filtered;
};

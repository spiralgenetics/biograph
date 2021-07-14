
#pragma once

#include "modules/bio_base/struct_var.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/io/transfer_object.h"
#include "modules/mapred/mapper.h"

struct filter_assemblies_params 
{
	TRANSFER_OBJECT 
	{
		VERSION(0); 
		FIELD(kmer_db);   // For kmer alignment
	} 

	std::string kmer_db;

	void validate() {}
}; 

class filter_assemblies_mapper: public typed_mapper<filter_assemblies_mapper,
	seq_position, struct_var, 
	seq_position, struct_var>
{
public:
	filter_assemblies_mapper(const std::string& params);
	~filter_assemblies_mapper();

	void set_watchdog(const std::function<void()>& watchdog) override
	{ 
		m_watchdog = watchdog; 
	}

	void setup() override;
	virtual void typed_map(const seq_position& key, const struct_var& value);	
	void install_metadata(meta::data& metadata) override;

private:
	std::function<void()> m_watchdog;
	filter_assemblies_params m_params;
	std::unique_ptr<const kmer_set> m_kdb;
	
	unsigned long m_mapped_count = 0;
	unsigned long m_filtered_count = 0;
	std::vector<unsigned> m_mismatch_distribution;
};

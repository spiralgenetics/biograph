
#pragma once

#include "modules/mapred/mapper.h"
#include "modules/bio_base/corrected_read.h"

class expand_mapper : public typed_mapper<expand_mapper, 
	std::string, corrected_read, 
	dna_sequence, uint64_t>
{
public:
	expand_mapper(const std::string& params) {}
	~expand_mapper() {}
	
	void typed_map(const std::string& key, const corrected_read& cr);	
};



#pragma once

#include "modules/mapred/mapper.h"
#include "modules/bio_base/dna_sequence.h"

class rotate_mapper : public typed_mapper<rotate_mapper, 
	dna_sequence, uint64_t, dna_sequence, int>
{
public:
	rotate_mapper(const std::string& params) {}
	~rotate_mapper() {}
	
	void typed_map(const dna_sequence& key, const uint64_t& val);	
};


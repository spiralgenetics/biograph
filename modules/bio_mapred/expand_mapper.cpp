
#include "modules/bio_mapred/expand_mapper.h"

REGISTER_1(mapper, expand, const std::string&);

void expand_mapper::typed_map(const std::string& key, const corrected_read& cr) 
{
	size_t size = cr.corrected.size();
	for(size_t i = 0; i < size; i++)
		output(cr.corrected.subseq(i, size - i), 1);
	dna_sequence rev = cr.corrected.rev_comp();
	for(size_t i = 0; i < size; i++)
		output(rev.subseq(i, size - i), 1);
}


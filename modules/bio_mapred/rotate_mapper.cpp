
#include "modules/bio_mapred/rotate_mapper.h"
#include "modules/io/log.h"

REGISTER_1(mapper, rotate, const std::string&);

void rotate_mapper::typed_map(const dna_sequence& key, const uint64_t& value) 
{
	//SPLOG("DOING ROTATE: %s, %d", key.as_string().c_str(), (int) value);
	output(key, -1);
	if (key.size() > 1)
		output(key.subseq(1, key.size() - 1), (int) key[0]);
}


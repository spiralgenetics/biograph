
#include "modules/bio_mapred/dna_sorter.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/io/msgpack_transfer.h"

REGISTER_1(sorter, dna, const std::string&);

int dna_sorter::compare(const std::string& key1, const std::string& key2) const
{
	dna_sequence s1;
	dna_sequence s2;
	if (key1 != "") msgpack_deserialize(s1, key1);
	if (key2 != "") msgpack_deserialize(s2, key2);

	if (s1 < s2) return -2;
	if (s2 < s1) return 2;
	return 0;
}

size_t dna_sorter::partition(const std::string& key, size_t num_partitions) const
{
	size_t tot = 0;
	if (num_partitions == 1) return 0;
	dna_sequence s1;
	msgpack_deserialize(s1, key);
	for(size_t i = 0; i < s1.size(); i++)
	{
		tot *= 5;
		tot += (int) s1[i];
	}
	return tot % num_partitions;
}


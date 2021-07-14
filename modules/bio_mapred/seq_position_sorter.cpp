#include "modules/io/msgpack_transfer.h"
#include "modules/bio_base/seq_position.h"
#include "modules/bio_mapred/seq_position_sorter.h"

REGISTER_1(sorter, seq_position, const std::string&);

int seq_position_sorter::compare(const std::string& key1, const std::string& key2) const
{
	seq_position s1;
	seq_position s2;
	msgpack_deserialize(s1, key1);
	msgpack_deserialize(s2, key2);

	if (s1.scaffold_id < s2.scaffold_id) return -2;
	if (s2.scaffold_id < s1.scaffold_id) return 2;
	if (s1.position < s2.position) return -1;
	if (s2.position < s1.position) return 1;
	return 0;
}

size_t seq_position_sorter::partition(const std::string& key, size_t num_partitions) const
{
	if (num_partitions == 1) return 0;
	seq_position s1;
	msgpack_deserialize(s1, key);
	return s1.scaffold_id % num_partitions;
}

#include "modules/mapred/pair_sorter.h"
#include "modules/bio_base/unaligned_read.h"

REGISTER_1(sorter, pair, const std::string&);

int pair_sorter::compare(const std::string& serialized_key1, const std::string& serialized_key2) const
{
	read_id	key1;
	read_id key2;
	msgpack_deserialize(key1, serialized_key1);
	msgpack_deserialize(key2, serialized_key2);
	std::string id1{key1.pair_name};
	std::string id2{key2.pair_name};
	
	if (id1 == id2) {
		return 0;
	}
	
	if (id1.empty() || id2.empty()) {
		return id1.empty() ? -2 : 2;
	}
	
	if (*(id1.crbegin()) != '1' && *(id1.crbegin()) != '2') {
		return id1 < id2 ? -2 : 2;
	}
	
	if (*(id2.crbegin()) != '1' && *(id2.crbegin()) != '2') {
		return id1 < id2 ? -2 : 2;
	}
	if (id1.substr(0, id1.size() - 1) == id2.substr(0, id2.size() - 1)) {
		return *(id1.crbegin()) < *(id2.crbegin()) ? -1 : 1;
	}
	
	return id1.substr(0, id1.size() - 1) < id2.substr(0, id2.size() - 1) ? -2 : 2;
}


size_t pair_sorter::partition(const std::string& key, size_t num_partitions) const
{
	size_t tot = 0;
	if (num_partitions == 1) {
		return 0;
	}
	
	if (key.empty()) {
		return 0;
	}

	for(size_t i = 0; i < key.size() - 1; i++)
	{
		tot *= 53;
		tot += key[i];
	}
	return tot % num_partitions;
}

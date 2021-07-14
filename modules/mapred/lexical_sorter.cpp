
#include "modules/mapred/lexical_sorter.h"

REGISTER_1(sorter, lexical, const std::string&);

int lexical_sorter::compare(const std::string& key1, const std::string& key2) const
{
	int comp = key1.compare(key2);
	return comp == 0 ? 0 : (comp < 0 ? -2 : 2);
}

size_t lexical_sorter::partition(const std::string& key, size_t num_partitions) const
{
	size_t tot = 0;
	if (num_partitions == 1)
		return 0;

	for(size_t i = 0; i < key.size(); i++)
	{
		tot *= 53;
		tot += key[i];
	}
	return tot % num_partitions;
}



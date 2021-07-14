
#ifndef __lexical_sorter_h__
#define __lexical_sorter_h__

#include "modules/mapred/sorter.h"

class lexical_sorter : public sorter
{
public:
	lexical_sorter(const std::string& params) {}
	int compare(const std::string& key1, const std::string& key2) const override;
	size_t partition(const std::string& key, size_t num_partitions) const override;
};

#endif


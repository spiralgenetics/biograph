#pragma once

#include "modules/mapred/sorter.h"

class pair_sorter : public sorter
{
public:
	pair_sorter(const std::string& params) {}

	int compare(const std::string& key1, const std::string& key2) const override;
	size_t partition(const std::string& key, size_t num_partitions) const override;
};

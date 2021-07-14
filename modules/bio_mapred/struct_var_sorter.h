#pragma once

#include "modules/mapred/sorter.h"

class struct_var_sorter : public sorter
{
public:
        struct_var_sorter(const std::string& params) {}
        int compare(const std::string& key1, const std::string& key2) const override;
        size_t partition(const std::string& key, size_t num_partitions) const override;
};

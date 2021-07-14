
#ifndef __dna_sorter_h__
#define __dna_sorter_h__

#include "modules/mapred/sorter.h"

class dna_sorter : public sorter
{
public:
        dna_sorter(const std::string& params) {}
        int compare(const std::string& key1, const std::string& key2) const override;
        size_t partition(const std::string& key, size_t num_partitions) const override;
};

#endif


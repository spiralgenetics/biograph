
#include "modules/mapred/sorter.h"

DEFINE_REGISTRY_1(sorter, std::string const&);

typedef simple_sorter<uint64_t> uint64_sorter;

REGISTER_1(sorter, uint64, const std::string&);

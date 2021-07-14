#include "modules/mapred/splitter.h"
#include "modules/io/registry.h"

DEFINE_REGISTRY_1(splitter, std::string const&);

REGISTER_1(splitter, null, const std::string&);

// Out of lined destructor to avoid missing vtable link error.
splitter::~splitter()
{
}

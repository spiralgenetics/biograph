
#include "modules/mapred/value_count_mapper.h"

void value_count_mapper::map(const std::string& key, const std::string& value, kv_sink& context)
{
	uint64_t one = 1;
	context.write(value, msgpack_serialize(one));	
}

REGISTER_1(mapper, value_count, const std::string&);



#ifndef __identity_mapper_h__
#define __identity_mapper_h__

#include "modules/mapred/mapper.h"

class identity_mapper : public mapper
{
public:
	identity_mapper(const std::string& params) {}
	void map(const std::string& key, const std::string& value, kv_sink& context) override
	{ context.write(key, value); }
};

#endif


#ifndef __read_qual_mapper_h__
#define __read_qual_mapper_h__

#include "modules/mapred/mapper.h"

class read_qual_mapper : public mapper
{
public:
	read_qual_mapper(const std::string& params) {}
        ~read_qual_mapper() {}
        
        void map(const std::string& key, const std::string& value, kv_sink& context) override;
};

#endif


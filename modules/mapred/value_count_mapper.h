
#ifndef __value_count_mapper_h__
#define __value_count_mapper_h__

#include "modules/mapred/mapper.h"

class value_count_mapper : public mapper
{
public:
	value_count_mapper(const std::string& params) {}
  void map(const std::string& key, const std::string& value, kv_sink& context) override;
};

#endif

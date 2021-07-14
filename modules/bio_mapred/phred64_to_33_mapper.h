#pragma once

#include "modules/mapred/mapper.h"

struct unaligned_read;
class phred64_to_33_mapper : public mapper
{
public:
	phred64_to_33_mapper(const std::string& params) {}
	virtual ~phred64_to_33_mapper() {}
        
	void map(const std::string& key, const std::string& value, kv_sink& context) override;

private:
	void convert_64_to_33(std::string& qualities_string) const;
	std::string make_exception_string(const unaligned_read& a_read, char bad_quality) const;
};

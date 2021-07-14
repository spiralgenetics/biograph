#pragma once

#include "modules/io/keyvalue.h"
#include "modules/mapred/path.h"

class file_info;

class file_info_reader : public kv_source
{
public:
	file_info_reader(const file_info& fi, const std::string& encoding);

	inline 
	std::string get_first_key()
	{
		return m_first_key;
	}

	bool read(std::string& key, std::string& value) override;
	
private:
	std::string                m_first_key;
	path                       m_path;
	std::string                m_encoding;
	std::unique_ptr<kv_reader> m_reader;
	std::unique_ptr<readable>  m_raw;
	std::unique_ptr<readable>  m_decoded;
};

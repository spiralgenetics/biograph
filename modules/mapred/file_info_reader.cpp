#include "modules/mapred/file_info_reader.h"
#include "modules/mapred/manifest.h"
#include "modules/io/keyvalue.h"
#include "modules/io/make_unique.h"

file_info_reader::file_info_reader(const file_info& fi, const std::string& encoding)
	: m_first_key(fi.first_key)
	, m_path(fi.file)
	, m_encoding(encoding)
{
	//SPLOG("in file_info_reader: using encoding: %s", encoding.c_str());
} 

bool file_info_reader::read(std::string& key, std::string& value)
{
	if (!m_reader) {
		m_raw = m_path.read();
		m_decoded = make_decoder(m_encoding, *m_raw);
		m_reader = make_unique<kv_reader>(*m_decoded);
	}
	return m_reader->read(key, value);
}

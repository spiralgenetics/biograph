
#include "modules/io/inverse_kvread.h"

inverse_kvread::inverse_kvread(reset_kv_source& source)
	: m_source(source)
	, m_loop_write(m_loop)
{}

size_t inverse_kvread::read(char* buf, size_t len)
{
	while (m_loop.size() < len)  // Try to add more data
	{
		std::string key;
		std::string value;
		bool r = m_source.read(key, value);
		if (r)
			m_loop_write.write(key, value);
		else 
			break; // No more to add
	}
	return m_loop.read(buf, len);
}

void inverse_kvread::reset()
{
	m_source.reset();
	m_loop.clear();
}


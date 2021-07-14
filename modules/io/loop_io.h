
#ifndef __loop_io_h__
#define __loop_io_h__

#include "modules/io/io.h"
#include <deque>

class loop_io : public readable, public writable 
{
public:
	loop_io() {}
	size_t size() { return m_buffer.size(); }
	size_t read(char* buf, size_t len) override
	{
		if (len > size())
			len = size();
		for(size_t i = 0; i < len; i++)
		{
			buf[i] = m_buffer.front();
			m_buffer.pop_front();
		}
		return len;
	}

	void write(const char* buf, size_t len) override
	{
		for(size_t i = 0; i < len; i++)
			m_buffer.push_back(buf[i]);
	}
	void clear() { m_buffer.clear(); }

private:
	std::deque<char> m_buffer;
};

#endif


#include "modules/web/chunked_encoding.h"
#include "modules/io/utils.h" 
#include "modules/io/log.h" 
#include <string.h> 

const unsigned int BUFFER_SIZE = 64*1024;

chunked_encoding_writable::chunked_encoding_writable(writable& sink)
	: m_sink(sink)
{
	m_begin = new char[BUFFER_SIZE];
	m_current_pos = m_begin;
	m_end = m_begin + (BUFFER_SIZE - 1);
}

chunked_encoding_writable::~chunked_encoding_writable()
{
	if(m_begin)
	{
		delete[] m_begin;
		m_begin = nullptr;
		m_current_pos = nullptr;
		m_end = nullptr;
	}
}

void chunked_encoding_writable::write(const char* buf, size_t len)
{
	if( m_current_pos + len <= m_end)
	{
		memcpy( m_current_pos, buf, len);
		m_current_pos += len;
	}
	else
	{
		size_t buffer_content_length = m_current_pos - m_begin,
		       chunk_size = buffer_content_length + len;
		std::string chunk_size_str = printstring("%lX\r\n", chunk_size);
		m_sink.write( chunk_size_str.c_str(), chunk_size_str.size() );
		m_sink.write( m_begin, buffer_content_length );
		m_sink.write( buf, len );
		m_sink.write("\r\n", 2);
		m_current_pos = m_begin;
	}
}

void chunked_encoding_writable::flush()
{
	size_t buffer_content_length = m_current_pos - m_begin;
	std::string chunk_size_str = printstring("%lX\r\n", buffer_content_length);
	m_sink.write( chunk_size_str.c_str(), chunk_size_str.size() );
	m_sink.write( m_begin, buffer_content_length );
	m_sink.write("\r\n", 2);
	m_current_pos = m_begin;
}

void chunked_encoding_writable::close()
{
	flush();
	m_sink.write("0\r\n\r\n", 5);
	m_sink.close();
}

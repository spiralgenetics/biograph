#include "modules/io/bzip.h"

#include <cstring>

bzip_reader::bzip_reader(readable& source) 
	: m_source(source)
{
	std::memset(&m_stream, 0, sizeof(m_stream));
	auto retcode = BZ2_bzDecompressInit(&m_stream, 0, 0);
	if (retcode != BZ_OK) {
		throw io_exception("bzip_reader::start> BZ2_bzDecompressInit() failed");
	}
}

bzip_reader::~bzip_reader()
{
	BZ2_bzDecompressEnd(&m_stream);
}

int bzip_reader::base_read(char* buf, size_t len)
{
	if (m_eof) {
		return 0;
	}

	m_stream.next_out = buf;
	m_stream.avail_out = len;
	while (m_stream.avail_out > 0) {
		if (m_stream.avail_in == 0) {
			m_stream.next_in = m_buf.data();
			m_stream.avail_in = m_source.read(m_buf.data(), m_buf.size());
		}
		
		auto retcode = BZ2_bzDecompress(&m_stream);
		if (retcode != BZ_OK) {
			m_eof = true;
			if (retcode == BZ_STREAM_END) {
				return len - m_stream.avail_out;
			}
			throw io_exception("bzip_reader::base_read> BZ2_bzDecompress() failed");
		}
	}
	
	// Must have read it all to get here
	return len;	
}

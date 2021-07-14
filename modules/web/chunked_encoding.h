#ifndef __chunked_encoding_h__
#define __chunked_encoding_h__

#include "modules/io/io.h"

// writes its input data into the Transfer-Encoding: chunked format, as defined in HTTP v1.1
// Note that this implementation does not currently support trailing headers.
class chunked_encoding_writable : public writable
{
public:
	chunked_encoding_writable(writable& sink);
	~chunked_encoding_writable();

	// Writes to its provided sink are put in a BUFFER_SIZE-byte buffer before a new chunk is produced.
	// if len of buf is greater than remaining buffer space, a new chunk is emitted that contains
	// the internal buffer data followed by the content of buf.
	void write(const char* buf, size_t len) override;
	void flush() override;
	void close() override; // calls flush before closing. Closing emits the last chunk of size zero, per protocol.
private:
	writable& m_sink;
	char *m_begin, *m_current_pos, *m_end;
};

#endif //__chunked_encoding_h__

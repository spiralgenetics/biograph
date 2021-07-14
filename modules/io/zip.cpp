#include "modules/io/zip.h"
#include "modules/io/log.h"
#include "modules/io/utils.h"
#include "base/base.h"

#include <cstring>

namespace
{

#define CODE_CASE(code) case code: return #code;

const char* zerr(int code)
{
	switch (code) {
	CODE_CASE(Z_OK);
	CODE_CASE(Z_STREAM_END);
	CODE_CASE(Z_NEED_DICT);
	CODE_CASE(Z_ERRNO);
	CODE_CASE(Z_STREAM_ERROR);
	CODE_CASE(Z_DATA_ERROR);
	CODE_CASE(Z_MEM_ERROR);
	CODE_CASE(Z_BUF_ERROR);
	CODE_CASE(Z_VERSION_ERROR);
	default:
		return "Unknown";
	}
}

} // namespace (anonymous)

zip_reader::zip_reader(readable& source, progress_t& update)
    : m_source(source), m_tracker(update) {
  std::memset(&m_stream, 0, sizeof(m_stream));
  auto retcode = inflateInit2(&m_stream, 32 + 15);
  if (retcode != Z_OK) {
    throw io_exception(
        printstring("zlib_reader> inflateInit2() failed: %s", zerr(retcode)));
  }

  m_read_thread = std::async(std::launch::async, [this]() {
    try {
      run_read_thread();
    } catch(...) {
      // Make sure the main thread doesn't sit around waiting forever
      // if we get an exception.
      m_closing = true;
      m_read_buffer_avail.notify_all();
      throw;
    }
  });
}

int zip_reader::read_internal(char* buf, size_t len)
{
	if (m_eof) {
		// SPLOG("zip_reader::base_read> EOF");
		return 0;
	}
	// SPLOG("zip_reader::base_read> %zu bytes", len);

	m_stream.next_out = reinterpret_cast<Bytef*>(buf);
	m_stream.avail_out = len;

	// while we still have room to put the output (decompressed) data
	while (m_stream.avail_out > 0) {
		// if we don't have anything more to decompress
		if (m_stream.avail_in == 0) {
			// get more space to put input (compressed) data 
			m_stream.next_in = reinterpret_cast<Byte*>(m_buf.data());
			// read compressed data and put it in m_buf, tell zlib how much we added
			m_stream.avail_in = m_source.read(m_buf.data(), m_buf.size());
		} else if (m_stream.avail_in == 1) {
			// We need to refill the buffer when it drops to one byte
			// to avoid the bug that arises when we split a gzip header
			// over two chunks.  If we have a single byte and it's 31
			// when Z_STREAM_END is set we have no way of knowing whether
			// the next chunk has a header or we're at EOF, so we avoid
			// the problem by making sure we never have a single byte
			// remaining in the buffer.
			m_buf[0] = *m_stream.next_in;
			m_stream.next_in = reinterpret_cast<Byte*>(m_buf.data());
			m_stream.avail_in = m_source.read(m_buf.data() + 1, m_buf.size() - 1) + 1;
		}

		int retcode = inflate(&m_stream, Z_SYNC_FLUSH);
		
		// We read again in case we consumed all the data at EOF.  If that happens,
		// we would return with m_eof true, but that's not what we want to happen
		// if m_source has more data available.  If we are actually at EOF, then
		// the extra reads should be harmless.
		if (m_stream.avail_in == 0) {
			m_stream.next_in = reinterpret_cast<Byte*>(m_buf.data());
			m_stream.avail_in = m_source.read(m_buf.data(), m_buf.size());
		} else if (m_stream.avail_in == 1) {
			m_buf[0] = *m_stream.next_in;
			m_stream.next_in = reinterpret_cast<Byte*>(m_buf.data());
			m_stream.avail_in = m_source.read(m_buf.data() + 1, m_buf.size() - 1) + 1;
		}

		switch (retcode) {
		case Z_STREAM_END:
			// SPLOG("zip_reader::base_read> Z_STREAM_END avail_out: %d", m_stream.avail_out);
			// SPLOG("zip_reader::base_read> Z_STREAM_END avail_in:  %d", m_stream.avail_in);
			m_tracker.final_update(m_stream.total_in, m_stream.total_out);
			if (check_eof()) {
				m_eof = true;
				// SPLOG("zip_reader::base_read> Z_STREAM_END: %zu bytes processed", len - m_stream.avail_out);
				return len - m_stream.avail_out;
			}
			retcode = inflateReset(&m_stream);
			if (retcode != Z_OK) {
				throw io_exception(printstring(
					"zip_reader::base_read> inflateReset() failed: %s", zerr(retcode)
				));
			}
			continue;
		case Z_BUF_ERROR:
			SPLOG("zip_reader::base_read> Z_BUF_ERROR: %zu bytes processed", len - m_stream.avail_out);
			return len - m_stream.avail_out;
		case Z_OK:
			// SPLOG("zip_reader::base_read> Z_OK");
			break;
		default:
			throw io_exception(printstring(
				"zip_reader::base_read> inflate() failed: %s", zerr(retcode)
			));
		}
	}
	return len;
}

void zip_reader::run_read_thread() {
  // Double buffer decompresison.  We decompress into work_buffer in
  // the background, and then swap it with done_buffer.  Then clients
  // read from done_buffer while we decompress into the next
  // work_buffer.
  std::unique_ptr<char[]> work_buffer(new char[k_decompress_buf_size]);
  std::unique_ptr<char[]> done_buffer(new char[k_decompress_buf_size]);
  for (;;) {
    int n_read = read_internal(work_buffer.get(), k_decompress_buf_size);
    if (n_read < 0) {
      SPLOG("zlib_reader: read_internal returned %d", n_read);
      n_read = 0;
    }
    CHECK_LE(n_read, k_decompress_buf_size);

    std::unique_lock<std::mutex> l(m_mu);
    while (m_out_buffer) {
      if (m_closing) {
        return;
      }
      m_read_buffer_consumed.wait(l);
    }

    if (m_closing) {
      return;
    }

    std::swap(work_buffer, done_buffer);
    m_out_buffer = done_buffer.get();
    m_out_size = n_read;
    m_read_buffer_avail.notify_one();
    if (n_read == 0) {
      CHECK(m_eof);
      return;
    }
    l.unlock();

    m_tracker.update(m_stream.total_in, m_stream.total_out);
  }
}

zip_reader::~zip_reader()
{
  {
    std::lock_guard<std::mutex> l(m_mu);
    m_closing = true;
    m_read_buffer_consumed.notify_all();
  }
  close_read_thread();
	int retcode = inflateEnd(&m_stream);
	if (retcode != Z_OK) {
		SPLOG("zlib_reader::~zip_reader> inflateEnd() failed: %s", zerr(retcode));
	}
}

void zip_reader::close_read_thread() {
  if (m_read_thread.valid()) {
    // Propagate any exceptions that might have happened on the read
    // thread.
    m_read_thread.get();
  }
}

int zip_reader::base_read(char* buf, size_t len) {
  std::unique_lock<std::mutex> l(m_mu);
  // SPLOG("zip_reader::base_read> %zu bytes", len);

  while (!m_out_buffer) {
    if (m_closing) {
      l.unlock();
      close_read_thread();
      return 0;
    }
    m_read_buffer_avail.wait(l);
  }

  if (m_closing) {
    l.unlock();
    close_read_thread();
    return 0;
  }

  if (m_out_size == 0) {
    CHECK(m_eof);
    return 0;
  }

  l.unlock();

  size_t size_to_read = std::min(len, m_out_size);
  memcpy(buf, m_out_buffer, size_to_read);
  m_out_buffer += size_to_read;
  m_out_size -= size_to_read;

  if (m_out_size == 0) {
    l.lock();
    m_out_buffer = nullptr;
    m_read_buffer_consumed.notify_one();
  }

  return size_to_read;
}

bool zip_reader::check_eof()
{
	// This code comes from zlib, gzread.c, gz_look().
	// The idea is to check if the next two bytes look like a GZIP header
	if (m_stream.avail_in > 1 &&
		m_stream.next_in[0] == 31 && m_stream.next_in[1] == 139) {
		// we found a GZIP header, don't transition to EOF
		return false;
	}
	return true;
}

zip_writer::zip_writer(
	writable& sink,
	progress_t& update,
	int compression_level,
	int compression_strategy)
	: m_sink(sink)
	, m_tracker(update)
{
	std::memset(&m_stream, 0, sizeof(m_stream));
	auto retcode = deflateInit2(
		&m_stream,
		compression_level,
		Z_DEFLATED,
		15 + 16,  // 15 windowBits, 16 = 'add a gzip header'
		9,        // use maximum memory
		compression_strategy
	);
	if (retcode != Z_OK) {
		throw io_exception(printstring(
			"zlib_writer> deflateInit2() failed: %s", zerr(retcode)
		));
	}
}

zip_writer::~zip_writer()
{
	base_close();
}

int zip_writer::base_write(const char* buf, int len)
{
	m_stream.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(buf));
	m_stream.avail_in = len;
	int retcode = compress(Z_NO_FLUSH);
	if (retcode != Z_OK) {
		throw io_exception(printstring(
			"zlib_writer::base_writer> deflate() failed: %s", zerr(retcode)
		));
	}
	return len;
}

int zip_writer::base_close()
{
	// to make ::close idempotent.
	if (m_closed) {
		return 0;
	}

	int retcode = compress(Z_FINISH);
	if (retcode != Z_STREAM_END) {
		throw io_exception(printstring(
			"zlib_writer::base_close> deflate() failed: %s", zerr(retcode)
		));
	}

	retcode = deflateEnd(&m_stream);
	if (retcode != Z_OK) {
		throw io_exception(printstring(
			"zlib_writer::base_close> deflateEnd() failed: %s", zerr(retcode)
		));
	}

	m_sink.close();
	m_closed = true;
	return 0;
}

int zip_writer::compress(int flush)
{
	int retcode = Z_OK;
	while ((m_stream.avail_in || flush == Z_FINISH) && retcode == Z_OK) {
		m_stream.next_out = reinterpret_cast<Bytef*>(m_buf.data());
		m_stream.avail_out = m_buf.size();

		retcode = deflate(&m_stream, flush);
		m_sink.write(m_buf.data(), m_buf.size() - m_stream.avail_out);
		m_tracker.update(m_stream.total_in, m_stream.total_out);
	}
	return retcode;
}

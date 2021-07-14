#include "base/base.h"
#include "modules/io/io.h"
#include "modules/io/utils.h"

#include <future>
#include <stdarg.h>
#include <memory.h>
#include <malloc.h>
#include <boost/lexical_cast.hpp>

const static size_t k_bufsize = 64 * 1024;

void writable::print(const char* format, ...)
{
	va_list vl;
	va_start(vl, format);

	char* buf;
	int len = vasprintf(&buf, format, vl);
	write(buf, len);
	free(buf);
	va_end(vl);
}

void writable::print(const std::string& s) {
  write(s);
}

void writable::write(const std::string& s) {
	write(s.c_str(), s.size());
}

bool readable::readline(std::string& line, size_t maxsize)
{
	char c;
	line.clear();
	while(line.size() < maxsize)
	{
		size_t r = this->read(&c, 1);
		if (r == 0) {
			if (line.size() == 0) {
				return false;
			}
			else {
				// last line must not have had a trailing \n
				return true;
			}
		}
		if (c == '\r')
			continue;
		if (c == '\n')
			return true;
		line.push_back(c);
	}
	throw io_exception("Line overflow in readable::readline");
}

bool readable::readline_no_copy(const char*& line, size_t& len, size_t maxlen) {
  bool result = readline(m_line_buffer, maxlen);
  line = m_line_buffer.data();
  len = m_line_buffer.size();
  return result;
}

read_wrapper::read_wrapper()
	: m_buf(new char[k_bufsize])
	, m_start(0)
	, m_end(0)
{}

read_wrapper::~read_wrapper()
{
	delete[] m_buf;
}

size_t read_wrapper::read(char* buf, size_t len)
{
	size_t tot_read = 0;

	/* First use pushback buffer if any */
	size_t read_from_buf = m_end - m_start;
	if (read_from_buf) {
		if (read_from_buf > len) {
			read_from_buf = len;
		}
		memcpy(buf, m_buf + m_start, read_from_buf);
		len -= read_from_buf;
		buf += read_from_buf;
		tot_read += read_from_buf;
		m_start += read_from_buf;
	}

	while (len) {
		const size_t max = static_cast<size_t>(std::numeric_limits<int>::max());
		ssize_t r = base_read(buf, std::min(len, max));
		if (r < 0) {
			throw io_exception(printstring("IO error on read: %ld, %d", r, errno));
		}
		if (r == 0) {
			return tot_read;
		}
		tot_read += r;
		buf += r;
		len -= r;
	}

	return tot_read;
}

bool read_wrapper::readline(std::string& line, size_t maxlen)
{
	line.clear();
	while(line.size() < maxlen)
	{
		if (m_start == m_end)
		{
			int r = base_read(m_buf, k_bufsize);
			if (r < 0)
			{
				throw io_exception("IO error on readline");
			}
			if (r == 0)
			{
				if (line.size() == 0) {
					return false;
				}
				else {
					// last line must not have had a trailing \n
					return true;
				}
			}
			m_start = 0;
			m_end = r;
		}
		// Position in m_buf to copy from once we find a line.
		size_t copy_start_pos = m_start;
		// Max position we can copy to before overrunning the maximum line size.
		size_t max_len_pos = m_start + (maxlen - line.size());
		// Max position we can copy to before we overrun maximum line size or
		// run out of buffer to copy.
		size_t max_end_pos = std::min(max_len_pos, m_end);
		CHECK_GT(max_end_pos, 0);

		while (m_start < max_end_pos) {
			char next = m_buf[m_start];
			if (next == '\n') {
				line.append(m_buf + copy_start_pos, m_start - copy_start_pos);
				// Advance past new line.
				m_start++;
				// End of line; return it.
				return true;
			}
			if (next == '\r') {
				// Found carriage return; ignore it.
				// Copy our line so far
				line.append(m_buf + copy_start_pos, m_start - copy_start_pos);
				// Skip \r
				m_start++;
				// Reset where we're copying from to
				// the character after the \r
				copy_start_pos = m_start;
				// And continue searching for a newline.
				continue;
			}
			// Continue searching for the end of the block to copy.
			m_start++;
		}
		line.append(m_buf + copy_start_pos, m_start - copy_start_pos);
	}
	throw io_exception(printstring("line too long: %lu bytes [max is %lu bytes]", line.size(), maxlen));
	// we used to have this:
	// throw io_exception("Long line ("+std::to_string(line.size())+" bytes) in readline.  maxlen is "+std::to_string(maxlen)+" bytes");
	// But a bug in MinGW 4.6.1 for to_string prevent that code from compiling
	// Bug should be fixed in 4.8.x
}

bool read_wrapper::readline_no_copy(const char*& line, size_t& len,
                                    size_t maxlen) {
  for (;;) {
    CHECK_GE(m_end, m_start);
    char* nl = (char*)memchr(m_buf + m_start, '\n', m_end - m_start);
    if (nl) {
      line = m_buf + m_start;
      len = nl - line;

      m_start = nl + 1 - m_buf;

      if (len && line[len - 1] == '\r') {
        len--;
      }
      if (len >= maxlen) {
        throw io_exception(
            printstring("line of length %lu exceeds %lu bytes in length", len, maxlen));
      }
      return true;
    }

    if ((m_end - m_start) >= maxlen) {
      throw io_exception(printstring(
          "line exceeds %lu bytes in length without a newline", maxlen));
    }

    memmove(m_buf, m_buf + m_start, m_end - m_start);
    m_end -= m_start;
    m_start = 0;

    if (m_end >= k_bufsize) {
      throw io_exception(
          printstring("Read buffer of %lu too small in readline", maxlen));
    }

    int nread = base_read(m_buf + m_end, k_bufsize - m_end);
    if (nread < 0) {
      throw io_exception("IO error on readline");
    }
    if (nread == 0) {
      len = m_end - m_start;
      line = m_buf + m_start;
      if (m_start == m_end) {
        return false;
      }
      return true;
    }

    m_end += nread;
    CHECK_LE(m_end, k_bufsize);
  }
}

void write_wrapper::write(const char* buf, size_t len)
{
	while(len > 0)
	{
		ssize_t r = base_write(buf, std::min(len, static_cast<size_t>(std::numeric_limits<int>::max())));
		if (r <= 0)
		{
			std::string error_string("IO error during write: ");
			error_string += std::to_string(r);
			throw io_exception(error_string);
		}
		len -= r;
		buf += r;
	}
}

void write_wrapper::close()
{
	int ret = base_close();
	if (ret < 0)
	{
		throw io_exception(printstring("Error (%d) in write_wrapper::close", ret));
	}
}

void write_wrapper::flush()
{
	int ret = base_flush();
	if (ret < 0)
	{
		throw io_exception(printstring("Error (%d) in write_wrapper::flush", ret));
	}
}

void io_copy(readable& source, writable& sink)
{
	std::array<char, 64*1024> buf;
	while (true) {
		size_t len = source.read(buf.data(), buf.size());
		if (len == 0) {
			break;
		}
		sink.write(buf.data(), len);
	}
}

void io_copy(readable& source, writable& sink, size_t max)
{
	size_t total = 0;
	std::array<char, 64*1024> buf;
	while (total < max) {
		size_t limit = std::min(max - total, buf.size());
		size_t len = source.read(buf.data(), limit);
		if (len == 0) {
			break;
		}
		sink.write(buf.data(), len);
		total += len;
	}
}

void io_copy_pairs(std::vector<io_pair_t> pairs)
{
	std::vector<std::future<void>> futures;

	for (const auto& pair : pairs) {
		futures.emplace_back(std::async(std::launch::async, [&] {
			io_copy(*pair.first, *pair.second);
		}));
	}

	for (auto& future : futures) {
		future.get();
	}
}

bool io_match(readable& in1, readable& in2, size_t& first_diff_pos, io_match_update_t& update)
{
	size_t len1, len2, modulo(1);
	char byte1, byte2;

	first_diff_pos = 0;
	while(1)
	{
		len1 = in1.read(&byte1, 1);
		len2 = in2.read(&byte2, 1);

		if (first_diff_pos % modulo == 0)
			modulo = update(first_diff_pos);

		if ( len1 == 0 && len2 == 0 ) break;
		if ( len1 != 0 && len2 == 0 ) return false;
		if ( len2 != 0 && len1 == 0 ) return false;

		if ( byte1 != byte2 ) return false;
		first_diff_pos++;
	}
	return true;
}

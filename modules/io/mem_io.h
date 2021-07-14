#pragma once

#include "base/base.h"
#include "modules/io/io.h"
#include "modules/io/track_mem.h"

#include <vector>
#include <algorithm>
#include <string.h>

// First, you write stuff into it.
// Then, you read from it, up to however many bytes to wrote in.
//
// It allows you to by-pass the filesystem when testing small data readers/writers.
// In particular, its content can be checked using test/compare_readables.h

class mem_io : public reset_readable, public writable 
{
public:
  mem_io() = delete;
	virtual ~mem_io() {
		if (m_buffer) {
          m_alloc.deallocate(m_buffer, m_allocated);
        }
	}
        mem_io(mem_io&& old) :
            m_alloc(old.m_alloc),
		m_offset(old.m_offset),
		m_buffer(old.m_buffer),
		m_allocated(old.m_allocated),
		m_size(old.m_size) {
		// Make sure the old buffer doesn't get freed again.
		old.m_buffer = nullptr;
        old.m_allocated = 0;
	}
  mem_io(const mem_io& old) : m_alloc(old.m_alloc) {
		write(old.m_buffer, old.m_size);
		m_offset = old.m_offset;
	}
	mem_io& operator=(const mem_io& old) {
		clear();
		write(old.m_buffer, old.m_size);
		m_offset = old.m_offset;
		return *this;
	}

  explicit mem_io(const std::string& str, const track_alloc& alloc) : m_alloc(alloc)
	{
		write(str.data(), str.size());
	}

	size_t read(char* buf, size_t len) override
	{
		len = std::min(len, m_size - m_offset);
		memcpy(buf, m_buffer + m_offset, len);
		m_offset += len;
		return len;
	}

	void write(const char* buf, size_t len) override
	{
		reserve(m_size + len);
		memcpy(m_buffer + m_size, buf, len);
		m_size += len;
	}

	size_t size() const
	{ 
		return m_size;
	}

	void reset() override
	{ 
		m_offset = 0;
	}

	void clear() 
	{ 
		m_size = 0;
		reset();
	}

	std::string str() const
	{
		return std::string(m_buffer, m_size);
	}

	void reserve(size_t size) {
		if (size <= m_allocated) {
			return;
		}
        size_t orig_allocated = m_allocated;
		while (m_allocated < size) {
			if (!m_allocated) {
				m_allocated = 1024;
			} else {
				m_allocated *= 2;
			}
		}
        char* old_buffer = m_buffer;
        m_buffer = m_alloc.allocate(m_allocated);
        if (old_buffer) {
          memcpy(m_buffer, old_buffer, orig_allocated);
          m_alloc.deallocate(old_buffer, orig_allocated);
        }
		CHECK(m_buffer);
	}
	void resize(size_t size) {
		reserve(size);
		m_size = size;
	}
	
	char* buffer() { return m_buffer; }
	const char* buffer() const { return m_buffer; }
	
private:
    track_mem::allocator<char> m_alloc;
	size_t m_offset = 0;
	char* m_buffer = nullptr;
	size_t m_allocated = 0;
	size_t m_size = 0;
};

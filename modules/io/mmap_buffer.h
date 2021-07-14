#pragma once

#include <string>

#include "base/base.h"
#include "modules/io/file_wrapper.h"
#include "modules/io/membuf.h"

class mmap_buffer : public mutable_membuf_impl
{
public:
	enum class mode { read_only, read_write, copy_on_write, read_populate }; 

	// Empty buffer
	mmap_buffer() = default;

	mmap_buffer(const mmap_buffer& rhs) = delete;
	mmap_buffer& operator=(const mmap_buffer&) = delete;
	
	mmap_buffer(mmap_buffer&& move_me) noexcept
	{
		std::swap(m_file, move_me.m_file);
		std::swap(m_buffer, move_me.m_buffer);
		std::swap(m_size, move_me.m_size);
		std::swap(m_uuid, move_me.m_uuid);
	}
	
	mmap_buffer& operator=(mmap_buffer&& rhs) noexcept
	{
		std::swap(m_file, rhs.m_file);
		std::swap(m_buffer, rhs.m_buffer);
		std::swap(m_size, rhs.m_size);
		std::swap(m_uuid, rhs.m_uuid);
		
		return *this;
	}

	// Create, always read_write
	mmap_buffer(const std::string& path, size_t size); 

	// Open existing
	mmap_buffer(const std::string& path, mode the_mode = mode::read_only);  

	~mmap_buffer();

	bool is_open() const; 

	// Create, always read_write
	void open(const std::string& path, size_t size); 

	// Open existing
	void open(const std::string& path, mode the_mode = mode::read_only); 

	void sync();
	void close();

	size_t size() override { return m_size; }
	const char* data() override { return m_buffer; }
	char* mutable_data() override {
	  CHECK(m_mode != mode::read_only);
	  CHECK(m_mode != mode::read_populate);
	  return m_buffer;
	}
	char* buffer() { return m_buffer; }
	const char* buffer() const { return m_buffer; }
	size_t size() const { return m_size; }
	const std::string& path() const { return m_file.path(); }

	std::string get_uuid() const { return m_uuid; }
	void set_uuid(const std::string& uuid) { m_uuid = uuid; }

  void truncate(size_t size);

private:
	mode m_mode;
	file_wrapper m_file;
	char* m_buffer = nullptr;
	size_t m_size = 0;
	std::string m_uuid;
};

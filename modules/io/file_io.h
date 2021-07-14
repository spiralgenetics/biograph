#pragma once

#include "modules/io/io.h"
#include <sys/types.h>
#include <stdio.h>

class file_reader : public read_wrapper
{
public:
	file_reader();
	file_reader(FILE* file);
	file_reader(const std::string& filename);
	file_reader(const file_reader& rhs) = delete;
	file_reader(file_reader&& rhs);
	file_reader& operator=(const file_reader& rhs) = delete;
	file_reader& operator=(file_reader&& rhs);
	~file_reader();
	void seek(uint64_t off);
	uint64_t size() const;
	size_t pos() const;
	void close() override;

private:
	int base_read(char* buf, size_t len) override;

	bool m_user_file;
	FILE* m_file;
	mutable int64_t m_size = -1; // Cache file size
};

class file_writer : public write_wrapper
{
public:
	file_writer(FILE* file);
	file_writer(const std::string& filename, bool append = false);
	file_writer(const file_reader& rhs) = delete;
	file_writer(file_writer&& rhs);
	file_writer& operator=(const file_writer& rhs) = delete;
	file_writer& operator=(file_writer&& rhs);
	~file_writer();
	void seek(uint64_t off);
	size_t pos() const;

private:
	int base_write(const char* buf, int len) override;
	int base_flush() override;
	int base_close() override;

	bool m_user_file;
	FILE* m_file;
	std::string m_filename;
};

// given a filename, it returns its size! It's magical!
off_t fsize(const std::string& filepath);

// given a filename, get the file contents as std::string
std::string slurp_file(const std::string& filepath);

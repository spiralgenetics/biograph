#pragma once

#include "modules/io/io.h"

#include <array>

extern "C"
{
#include <bzlib.h>
};

class bzip_reader : public read_wrapper
{
public:
	bzip_reader(readable& source);
	~bzip_reader();

private:
	int base_read(char* buf, size_t len) override;

	readable& m_source;
	bool m_eof = false;

	bz_stream m_stream;
	std::array<char, 16*1024> m_buf;
};

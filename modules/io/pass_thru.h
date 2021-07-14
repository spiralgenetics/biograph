#pragma once

#include "modules/io/io.h" 
#include "modules/io/log.h" 

class pass_thru_writable : public writable
{
public:
	pass_thru_writable(writable& in)
		: m_in(in)
	{}

	void print(const char* format, ...) __attribute__((format(printf,2,3)));

	void write(const char* buf, size_t len) override
	{
		m_in.write(buf, len); 
	}

	void flush() override
	{
		m_in.flush(); 
	}

	void close() override
	{
		m_in.close(); 
	}

private:
	writable& m_in;
};

class pass_thru_readable : public readable
{
public:
	pass_thru_readable(readable& in)
		: m_in(in)	
	{}

	size_t read(char* buf, size_t len) override
	{
		return m_in.read(buf, len); 
	}

	bool readline(std::string& line, size_t maxlen) override
	{
		return m_in.readline(line, maxlen); 
	}

private:
	readable& m_in;
};

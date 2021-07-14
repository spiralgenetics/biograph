#pragma once

#include "modules/io/io.h"
#include "modules/io/utils.h"

#include <string.h>

class tab_writer
{
public:
	tab_writer(std::string& out) 
		: m_first(true)
		, m_out(out) 
	{}
	
	bool is_eof() { return false; }
	
	void write_tab()
	{
		if (m_first) {
			m_first = false;
		}
		else {
			m_out.push_back('\t');
		}
	}

	void write_unencoded(const char* str)
	{
		write_tab();
		for (size_t i = 0; i < strlen(str); i++) {
			m_out.push_back(str[i]);
		}
	}

	void write_char(char c)
	{
		if (c == '\\') { 
			m_out.push_back('\\'); 
			m_out.push_back('\\'); 
		}
		else if (c == '\t') { 
			m_out.push_back('\\'); 
			m_out.push_back('t'); 
		}
		else if (c == '\0') { 
			m_out.push_back('\\'); 
			m_out.push_back('0'); 
		}
		else {
			m_out.push_back(c);
		}
	}

	void write_number(const char prefix, const unsigned int& ir)
	{
		unsigned int i = ir;
		if (i == 0) {
			write_unencoded("0");
			return;
		}
		char buf[20];
		size_t pos = 20;
		buf[--pos] = 0;
		while (i != 0) {
			buf[--pos] = '0' + (i % 10);
			i /= 10;
		}
		if (prefix) {
			buf[--pos] = prefix;
		}

		write_unencoded(buf + pos);
	}
	void transfer(const unsigned int& i) { write_number(0, i); }
	void transfer(const int& i) { if (i < 0) write_number('-', -i); else write_number(0, i); }
	void transfer(const unsigned short& i) { write_number(0, i); }
	void transfer(const short& i) { if (i < 0) write_number('-', -i); else write_number(0, i); }
	void transfer(const unsigned char& i) { write_number(0, i); }
	void transfer(const char& c) 
	{ 
		write_tab();
		write_char(c);
	}
	void transfer(const bool& b) { write_number(0, (unsigned int) b); }
	void transfer(const std::string& str)
	{
		write_tab();
		for (size_t i = 0; i < str.size(); i++) {
			write_char(str[i]);
		}
	}
	template<class X>
	void transfer(const X& obj)
	{
		const_cast<X&>(obj).transfer(*this);
	}

private:
	bool m_first;
	std::string& m_out;
};

class tab_reader
{
public:
	tab_reader(const std::string& in) 
		: m_cur(0)
		, m_str(in) 
	{}

	char peek() const 
	{
		if (m_cur == m_str.size()) {
			return '\t'; 
		}
		return m_str[m_cur]; 
	}

	bool is_eof() const
	{
		return (m_cur == m_str.size());
	}

	char next() 
	{ 
		if (m_cur == m_str.size()) {
			throw io_exception("tab_reader: unexpected end of string"); 
		}
		return m_str[m_cur++];
	}

	void read_tab()
	{
		if (m_cur != 0 && next() != '\t') {
			throw io_exception("read_tab: expected tab");
		}
	}

	void skip_field()
	{
		read_tab();
		while (peek() != '\t') {
			next();
		}
	}

	int read_sign()
	{
		int sign;
		if (peek() == '-') {
			sign = -1;
			next();
		}   
		else {
			sign = 1;
		}

		return sign;
	}

	unsigned int read_value()
	{
		unsigned int value = 0;
		while (peek() != '\t') {
			value *= 10;
			char c = next();
			if (c < '0' || c > '9') {
				throw io_exception("read_signed_number: expected digit");
			}
			value += (c - '0');	
		}
		return value;
	}

	char read_char()
	{
		if (peek() == '\t' || peek() == 0) {
			throw io_exception("read_char: unexpected tab"); 
		}
		if (peek() == '\\') {
			next();
			char c = next();
			if (c == '\\') return '\\';
			if (c == 'N') return 0;  // Special placeholder for 'null' is SQL
			if (c == 't') return '\t';
			if (c == 'n') return '\n';
			if (c == '\t') return '\t';
			if (c == '0') return 0;
			throw io_exception(printstring("read_char: bad escape sequence: '%c'", c));
		}
		return next();
	}

	void transfer(size_t& i) { read_tab(); i = read_value(); }
	void transfer(unsigned int& i) { read_tab(); i = read_value(); }
	void transfer(int& i) { read_tab(); i = read_sign(); i *= read_value(); }
	void transfer(int& i, int def) 
	{
		read_tab();
		if (peek() == '\\') {
			next(); 
			char c = next();
			if (c == 'N') { 
				i = def; 
				return; 
			}
			throw io_exception("transfer: failed to find null");
		}
		i = read_sign(); 
		i *= read_value();
	}
	void transfer(unsigned short& i) { read_tab(); i = read_value(); }
	void transfer(short& i) { read_tab(); i = read_sign(); i *= read_value(); }
	void transfer(unsigned char& i) { read_tab(); i = read_value(); }
	void transfer(bool& i) { read_tab(); i = read_value(); }
	void transfer(char& i) 
	{ 
		read_tab();  
		i = read_char();
	}

	void transfer(std::string& s) 
	{
		read_tab();
		s.clear();
		while (peek() != '\t') {
			s.push_back(read_char());
		}
		if (s.size() == 1 && s[0] == 0) { // Handle 'null' strings by making them empty
			s.clear();
		}
	}

	void check_end()
	{	
		if (m_cur != m_str.size()) {
			throw io_exception("extra data in tabd");
		}
	}

	template<class X>
	void transfer(X& x)
	{
		x.transfer(*this);
	}
	
private:
	size_t m_cur;
	std::string m_str;
};

template<class Obj>
inline std::string tabd_serialize(const Obj& value)
{
	std::string out;
	tab_writer tw(out);
	tw.transfer(value);
	return out;
}

template<class Obj>
inline void tabd_deserialize(Obj& value, const std::string& str)
{
	tab_reader tr(str);
	tr.transfer(value);
	tr.check_end();
}


#define TABD_TRANSFER template<class T> void transfer(T& t) 
#define TFIELD(name) t.transfer(name)
#define TFIELD_DEF(name, def) t.transfer(name, def)
#define TFIELD_SKIP(name) t.skip_field()
#define TFIELD_MAYBE_SKIP(name) if (!t.is_eof()) t.skip_field()

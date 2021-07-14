
#include "modules/io/keyvalue.h"
#include <string.h>

size_t k_bufsize = 4096;

/* Helper functions that serialize/deserialize small, variable length 'size_t' elements 
 * each byte is 7 bits of the number, and the high bit is 0 if there are more bytes to follow, 1 if we are done
 */

// Returns false if at EOF, true otherwise, throws on err or EOF in middle of number
size_t read_small_size(readable& in)
{
	size_t r = 0;
	unsigned char next;
	do
	{
		if (in.read((char*) &next, 1) != 1)
			throw io_exception("Got end of file while reading size");
		r <<= 7;
		r |= (next & 0x7f);
	}
	while((next & 0x80)==0);
	return r;
}

void write_small_size(writable& out, size_t r)
{
	int c = 9;
	unsigned char buf[10];
	while(r >= 128)
	{
		buf[c--] = r & 0x7f;
		r >>= 7;
	}
	buf[c] = r;
	buf[9] |= 0x80;
	out.write((char *) (buf + c), 10 - c);
}
		

/* kv reader supports two formats: old and new
   Old format is 'K', key (no nulls), '\0', 'V', value (no nulls) '\0'
   New format is "2KVP" (header, only once), then per row: key small_size, key data, value small_size, value data
*/
kv_reader::kv_reader(readable& source) 
	: m_source(source)	
{}

kv_reader::~kv_reader()
{}

bool kv_reader::read(std::string& key, std::string& value)
{
	char first;
	if (m_source.read(&first, 1) == 0)
		return false;  // Empty file 
	if (first == 'N')
		read_newstyle(key, value);
	else if (first == 'K')
		read_oldstyle(key, value);
	else
		throw io_exception("Invalid kv-read start char (must be N or K)");

	return true;
}

void kv_reader::read_newstyle(std::string& key, std::string& value)
{
	size_t key_len, value_len;
	key_len = read_small_size(m_source);
	key.resize(key_len);
	if (m_source.read(&key[0], key_len) != key_len)
		throw io_exception("Unexpected EOF in key data");

	value_len = read_small_size(m_source);
	value.resize(value_len);
	if (m_source.read(&value[0], value_len) != value_len)
		throw io_exception("Unexpected EOF in value data");
}

void kv_reader::read_oldstyle(std::string& key, std::string& value)
{
	read_terminated(key);

	char magic;
	if (m_source.read(&magic, 1) != 1)
		throw io_exception("Key without value");
	if (magic != 'V')
		throw io_exception("Value missing magic id");
	read_terminated(value);
}

// TODO: This should probably be merged with readline in a generic fashion
void kv_reader::read_terminated(std::string& out)
{
	out.clear();
	while(true)
	{
		char c;
		if (m_source.read(&c, 1) != 1)
			throw io_exception("Unclean ending in read_terminated");
		if (c == 0)
			break;
		out.push_back(c);
	}
}

kv_writer::kv_writer(writable& sink)
	: m_sink(sink)
{}

kv_writer::~kv_writer() {}

void kv_writer::flush() {}
void kv_writer::close() { flush(); }

void kv_writer::write(const std::string& key, const std::string& value)
{
	m_sink.write("N", 1);
	write_small_size(m_sink, key.size());
	m_sink.write(&key[0], key.size());
	write_small_size(m_sink, value.size());
	m_sink.write(&value[0], value.size());
	/*  Old style
	m_sink.write("K", 1);
	m_sink.write(&key[0], key.size());
	m_sink.write("\0V", 2);
	m_sink.write(&value[0], value.size());
	m_sink.write("\0", 1);
	*/
}

static size_t compute_overhead(size_t sz)
{
	size_t r = 1;
        while(sz >= 128)
        {
                sz >>= 7;
		r++;
        }
	return r;
}

size_t kv_serial_size(size_t keysize, size_t valuesize)
{
	return 1 + compute_overhead(keysize) + compute_overhead(valuesize) + keysize + valuesize;
}


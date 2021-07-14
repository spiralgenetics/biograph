#pragma once

#include "modules/io/io.h"
#include "modules/io/msgpack_transfer.h"

class kv_source
{
public:
	virtual ~kv_source() = default;

	// Reads a key and value and returns true.  Returns false (and leaves
	// key and value alone) on EOF.  Throws io_exception for anything unexpected
	// TODO: should we allow user to specify max lengths?
	virtual bool read(std::string& key, std::string& value) = 0;

	// Nice version for 'typed' data
	template<class KeyType, class ValueType>
	bool read_msgpack(KeyType& key, ValueType& value)
	{
		std::string key_str, value_str;
		if (!read(key_str, value_str)) {
			return false;
		}
		msgpack_deserialize(key, key_str);
		msgpack_deserialize(value, value_str);
		return true;
	}
};

class reset_kv_source : public kv_source
{
public:
	virtual ~reset_kv_source() = default;
	virtual void reset() = 0;
};

class kv_sink
{
public:
	virtual ~kv_sink()
	{ 
		this->close(); 
	}

	// Write a key/value pair
	virtual void write(const std::string& key, const std::string& value) = 0;

	// Typed version
	template<class KeyType, class ValueType>
	void write_msgpack(const KeyType& key, const ValueType& value)
	{
		write(
			msgpack_serialize(key), 
			msgpack_serialize(value)
		);
	}
	
	virtual void close() {}
};

class kv_reader : public kv_source
{
public:
	kv_reader(readable& source);
	virtual ~kv_reader();

	bool read(std::string& key, std::string& value) override;

private:
	void read_oldstyle(std::string& key, std::string& value);
	void read_newstyle(std::string& key, std::string& value);
	void read_terminated(std::string& what);

	readable& m_source;
};

class kv_writer : public kv_sink
{
public:
	kv_writer(writable& sink);
	virtual ~kv_writer();

	void write(const std::string& key, const std::string& value) override;
	virtual void flush(); // write all buffered data to m_sink
	void close() override;

private:
	writable& m_sink;
};

size_t kv_serial_size(size_t keysize, size_t valuesize);

inline 
void kv_copy(kv_source& in, kv_sink& out)
{
	std::string key;
	std::string value;
	while (in.read(key, value)) {
		out.write(key, value);
	}
}


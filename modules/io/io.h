#pragma once

#include <boost/format.hpp>

#include <string>
#include <stdexcept>
#include <functional>
#include <vector>
#include <stdio.h>

class io_exception : public std::runtime_error
{
public:
	io_exception(const std::string& message) : std::runtime_error(message) {}
	io_exception(const boost::format& formatted_message) : std::runtime_error(formatted_message.str()) {}

	std::string message() const { return what(); }
};

class readable
{
public:
	virtual ~readable() {}

	// Always reads until len or EOF, returns total read, throws on err
	virtual size_t read(char* buf, size_t len) = 0;

	template<typename ToRead> size_t typed_read(ToRead* ptr, size_t count)
	{
		return read(reinterpret_cast<char*>(ptr), count * sizeof(ToRead));
	}

	// Reads till an EOF or EOL, return true if line was full (includes EOL), false if EOF first
	// In the case of return of false, line contains chars between last EOL and EOF
	// Never reads more than maxlen, throws if err or long line
	// Readable provides a (terrible) default implementation, read_wrapper's is better
	virtual bool readline(std::string& line, size_t maxlen);
	// Variant that returns a pointer into an internal buffer, so that
	// fewer copies are required.  This buffer may be overwritten on
	// the next call to readline() or close().
	virtual bool readline_no_copy(const char*& line, size_t& len, size_t maxlen);
	virtual void close() {}

 private:
  // Buffer for readline_no_copy to use for readables that don't
  // support no-copy readline.
  std::string m_line_buffer;
};

class reset_readable : public readable
{
public:
	virtual void reset() = 0;
};

/* A helper class to wrap objects with unix read semantics (-1 on err, not full reads, etc) */
class read_wrapper : public readable
{
public:
	read_wrapper();
	virtual ~read_wrapper();
	size_t read(char* buf, size_t len) override;
	void close() override {}
	bool readline(std::string& line, size_t maxlen) override;
	bool readline_no_copy(const char*& line, size_t& len, size_t maxlen) override;

private:
	virtual int base_read(char* buf, size_t len) = 0;
	char* m_buf;
	size_t m_start;
	size_t m_end;
};

class writable
{
public:
	virtual ~writable() {}

	// A fun helper function
	void print(const char* format, ...) __attribute__((format(printf,2,3)));

	// A useful helper function
	void print(const std::string& s);
	void write(const std::string& s);

	// Always does a full write, throws on errors, may do buffering prior to close
	virtual void write(const char* buf, size_t len) = 0;
	template<typename ToWrite> void typed_write(const ToWrite* ptr, size_t count)
	{
		write(reinterpret_cast<const char*>(ptr), count * sizeof(ToWrite));
	}

	// Flush data before returning
	virtual void flush() {}

	// Since you can write, it must make sense to 'finish' somehow, throws on error
	virtual void close() {}
};

/* A helper class to wrap objects with unix write semantics */
class write_wrapper : public writable
{
public:
	void write(const char* buf, size_t len) override;
	void flush() override;
	void close() override;

private:
	virtual int base_close() { return 0; }
	virtual int base_flush() { return 0; }
	virtual int base_write(const char *buf, int len) = 0;
};

// read 64KB at a time from in and writes it to out.
// does not close out
void io_copy(readable& source, writable& sink);
void io_copy(readable& source, writable& sink, size_t max);

// copy multiple pairs of io streams simultaneously
typedef std::pair<readable*, writable*> io_pair_t;
void io_copy_pairs(std::vector<io_pair_t> pairs);

typedef const std::function<size_t (size_t input_read)> io_match_update_t;
inline size_t io_match_no_update(size_t) { return (size_t)1; }

// Returns true if the two readables match byte for byte.
// Returns false as soon a difference is found:
// either two bytes differ in value,
// or one readable is exhausted before the other one.
//
// If false is returned, first_diff_pos will be set to the position
// of the first non-matching byte.
// That means, in1 and in2 matched byte-for-byte in the range [0..first_diff_pos-1].
//
// Implementation detail: the reading happens one byte at a time
// so it is slow.
// Which is fine because this function is only meant for testing purposes.
// If you can't tolerate its speed, please rewrite a buffered version.
//
bool io_match(
	readable& in1,
	readable& in2,
	size_t& first_diff_pos,
	io_match_update_t& update = io_match_no_update);

// A writable class that throws away any data written to it like /dev/null.
class null_writable : public writable
{
public:
	void write(const char* buf, size_t len) override {}
};

// multi_writer creates a writer that duplicates its writes to all the provided writers,
// similar to the Unix tee(1) command.
class multi_writer : public writable
{
public:
	template <class ...T>
	multi_writer(T... args)
		: m_writer({args...})
	{}

	void write(const char* buf, size_t len) override
	{
		for (const auto& writer : m_writer) {
			writer->write(buf, len);
		}
	}

	void flush() override
	{
		for (const auto& writer : m_writer) {
			writer->flush();
		}
	}

	void close() override
	{
		for (const auto& writer : m_writer) {
			writer->close();
		}
	}

private:
	std::vector<writable*> m_writer;
};

#pragma once

#include "modules/io/io.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/temp_file.h"
#include "modules/bio_base/dna_sequence.h"

class dna_sequence;

class dna_writer
{
public:
	explicit dna_writer(writable& target) : m_target(target) {}
	void write(const dna_sequence& seq_to_write);
	void flush() { m_target.flush(); }
	void close() { m_target.close(); }

private:
	writable& m_target;
	std::string m_packed_seq;
};

class dna_reader
{
public:
	dna_reader() = default;
	dna_reader(dna_reader&&) = default;
	dna_reader& operator=(dna_reader&&) = default;
	dna_reader(const dna_reader&) = delete;
	dna_reader& operator=(const dna_reader&) = delete;
	explicit dna_reader(std::unique_ptr<readable> source)
		: m_source(std::move(source))
	{}
	// Returns empty dna_sequence on EOF.
	dna_sequence read();

private:
	std::unique_ptr<readable> m_source;
	std::string m_packed_seq;
};

// Reads a DNA sequence in advance and buffers it for later use.
class dna_buffer
{
public:
	dna_buffer() = default;
	explicit dna_buffer(const std::string& file_path)
		: m_dna_reader(make_unique<file_reader>(file_path))
		, m_current_sequence(m_dna_reader.read())
	{}

	const dna_sequence& get_sequence() const { return m_current_sequence; }
	bool at_eof() const { return m_current_sequence.size() == 0; }
  virtual ~dna_buffer() = default;
	virtual void advance() { m_current_sequence = m_dna_reader.read(); }

protected:
	dna_reader m_dna_reader;
	dna_sequence m_current_sequence;
};

// Reads a DNA sequence sequentially from multiple files.  The sequences should not be
// broken across more than one file.  When a file is done, the reader moves on to the
// next file until files are exhausted.  The reader only advances, never retreats.
class multi_file_dna_buffer : public dna_buffer
{
	using temp_files_t = std::vector<std::shared_ptr<scoped_temp_file>>;
	
public:
	explicit multi_file_dna_buffer(temp_files_t temp_files)
		: dna_buffer()
		, m_source_files(std::move(temp_files))
		, m_current_file(m_source_files.begin())
	{ advance(); }
	
	void advance() override;
private:
	temp_files_t m_source_files;
	temp_files_t::const_iterator m_current_file;
};

#pragma once

#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_format/importer.h"
#include "modules/bio_format/exporter.h"

class fastq_reader : public kv_source
{
public:
	fastq_reader(readable& source, bool keep_quality = true);
	bool read(std::string& key, std::string& value) override;
	bool read(read_id& id, unaligned_read& value);
	bool read(read_id& id, unaligned_reads& value);

	size_t get_bases() const 
	{
		return m_bases;
	}

private:
	readable& m_source;
	int m_linenum = 0;
	size_t m_bases = 0;
  bool m_keep_quality = true;
};

class fastq_importer : public importer
{
public:
	fastq_importer(readable& source, bool /*unused*/, const std::string& /*unused*/)
		: m_reader(source) 
	{}

	fastq_importer(readable& source) 
		: m_reader(source) 
	{}
	
	void import(kv_sink& sink, simple_metadata& meta) override;


private:
	fastq_reader m_reader;
};

class writable;
class fastq_exporter : public exporter
{
public:
	fastq_exporter(writable& byte_sink) 
		: exporter(byte_sink) 
	{}

	fastq_exporter(writable& byte_sink, bool /*unused*/, const std::string& /*unused*/)
		: exporter(byte_sink) 
	{}

	void write(const std::string& key, const std::string& value) override;
	void write(const read_id& id, const unaligned_reads& reads);

};

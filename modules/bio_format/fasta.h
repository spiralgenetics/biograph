#pragma once

#include "modules/bio_format/importer.h"
#include "modules/bio_format/exporter.h"

// Default readline length. Some people have ridiculously long fasta lines;
// cut it off at 1GB.
static const int k_maxline = 1024*1024*1024;

class fasta_importer : public importer
{
public:
	fasta_importer(readable& source, bool /*unused*/, const std::string& /*unused*/)
		: m_source(source)
	{}
	fasta_importer(readable& source)
		: m_source(source)
	{}

	void import(kv_sink& sink, simple_metadata& meta) override;


private:
	readable& m_source;
};

class writable;
class fasta_exporter : public exporter
{
public:
	fasta_exporter(writable& byte_sink)
		: exporter(byte_sink)
	{}

	fasta_exporter(writable& byte_sink, bool /*unused*/, const std::string& /*unused*/)
		: exporter(byte_sink)
	{}

	void write(const std::string & key, const std::string & value) override;

};

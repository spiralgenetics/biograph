#pragma once

#include "modules/bio_format/exporter.h"

class kmer_count_exporter : public exporter
{
public:
	kmer_count_exporter(writable & byte_sink, unsigned kmer_size) 
		: exporter(byte_sink)
		, m_kmer_size(kmer_size) 
	{}

	kmer_count_exporter(writable & byte_sink, bool ignore, const std::string& params)
		: exporter(byte_sink)
		, m_kmer_size(stoi(params))
	{} 

	void write(const std::string& key, const std::string& value) override;
	
	
private:
	unsigned m_kmer_size;
};

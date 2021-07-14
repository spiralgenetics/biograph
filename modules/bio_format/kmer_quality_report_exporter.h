#pragma once

#include "modules/bio_format/exporter.h"

struct point
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(x, TF_STRICT); 
		FIELD(y, TF_STRICT); 
	}
	size_t x;
	size_t y; 
};

class kmer_quality_report_exporter : public exporter
{
public:
	kmer_quality_report_exporter(writable& byte_sink, const std::string& dataset_name);

	void write(const std::string& key, const std::string& value) override;
	void write_header() override;
	void write_footer() override;

private:
	std::string m_dataset_name;
	std::vector<point> m_data;
};

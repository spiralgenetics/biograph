#pragma once

#include "modules/bio_base/struct_var.h"
#include "modules/bio_format/exporter.h"

class corrected_reads_exporter : public exporter
{
public:
	corrected_reads_exporter(writable& byte_sink) 
		: exporter(byte_sink) 
	{}
	
	corrected_reads_exporter(writable& byte_sink,
		bool /*unused*/,
		const std::string& /*unused*/)
		: exporter(byte_sink)
	{}

	void write(const std::string& key, const std::string& value) override;
};


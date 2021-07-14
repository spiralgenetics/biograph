#pragma once

#include "modules/bio_base/struct_var.h"
#include "modules/bio_format/exporter.h"

class supporting_reads_exporter : public exporter
{
public:
	supporting_reads_exporter(writable& byte_sink) 
		: exporter(byte_sink) 
	{}
	
	supporting_reads_exporter(
		writable& byte_sink,
		bool /*unused*/,
		const std::string& /*unused*/)
		: exporter(byte_sink)
	{}

	void write(const std::string& key, const std::string& value) override;

};


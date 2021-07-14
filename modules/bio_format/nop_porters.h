#pragma once

#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_format/importer.h"
#include "modules/bio_format/exporter.h"
#include "modules/io/io.h"

class nop_importer : public importer
{
public:
	nop_importer(readable& source, bool /*unused*/, const std::string& /*unused*/)
		: m_source(source)
	{}

	void import(kv_sink& sink, simple_metadata& meta) override;

private:
	readable& m_source;
};

class nop_exporter : public exporter
{
public:
	nop_exporter(writable& byte_sink, bool /*unused*/, const std::string& /*unused*/)
		: exporter(byte_sink)
	{}

	void write(const std::string& key, const std::string& value) override;

};

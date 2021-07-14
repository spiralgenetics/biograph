#pragma once

#include "modules/io/keyvalue.h"
#include "modules/io/registry.h"

class exporter : public kv_sink
{
public:
	exporter(writable& byte_sink)	
		: m_sink(byte_sink)
	{}
	
	virtual ~exporter() = default;

	void write(const std::string& key, const std::string& value) override = 0;
	virtual void write_header() {}
	virtual void write_footer() {}

	void close() override
	{
		m_sink.close();
	}

	void export_from(kv_source& source)
	{
		write_header();

		std::string key;
		std::string value;
		while (source.read(key, value)) {   
			write(key, value);
		}   

		write_footer();
		close();
	}

protected:
	writable& m_sink;
};

DECLARE_REGISTRY_3(exporter, writable&, bool, std::string const&);

#define __STDC_FORMAT_MACROS

#include "modules/mapred/histogram_export.h"
#include "modules/io/msgpack_transfer.h"

#include <inttypes.h>

histogram_exporter::histogram_exporter(writable& byte_sink) 
	: exporter(byte_sink) 
{}

void histogram_exporter::write(const std::string& key, const std::string& value)
{
	uint64_t k, v;
	msgpack_deserialize(k, key);
	msgpack_deserialize(v, value);

	m_sink.print("%" PRIu64 "\t%" PRIu64 "\n", k, v);
};

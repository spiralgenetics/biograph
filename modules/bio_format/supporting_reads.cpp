#include "modules/bio_format/supporting_reads.h"
#include "modules/bio_base/struct_var.h"
#include "modules/io/registry.h"


REGISTER_3(exporter, supporting_reads, writable&, bool, const std::string&);

void supporting_reads_exporter::write(const std::string& key, const std::string& value)
{
	struct_var_key svkey;
	read_support support;
	msgpack_deserialize(svkey, key);
	msgpack_deserialize(support, value);
	
	std::string line = printstring("%u:%u\t%s\t%s\t%s\t%s\t%s\n",
		svkey.variation_id,
		svkey.read_id,
		support.name.c_str(),
		support.original.as_string().c_str(),
		support.corrected.as_string().c_str(),
		support.quality.c_str(),
		support.flipped ? "Flipped" : ""
	);
	
	m_sink.write(line.c_str(), line.size());
}

#include "modules/bio_format/nop_porters.h"
#include "modules/io/simple_metadata.h"
#include "modules/io/utils.h"
#include "modules/io/log.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/registry.h"


REGISTER_3(importer, nop, readable&, bool, const std::string&);
REGISTER_3(exporter, nop, writable&, bool, const std::string&);


void nop_importer::import(kv_sink& sink, simple_metadata& meta)
{
	SPLOG("nop_importer::import>");

	null_writable devnull;
	io_copy(m_source, devnull);

	SPLOG("nop_importer::import> flushed to /dev/null");
}

void nop_exporter::write(const std::string& key, const std::string& value)
{
}

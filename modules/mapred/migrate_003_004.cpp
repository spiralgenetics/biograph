#include "modules/mapred/migrate.h"
#include "modules/io/log.h" 

// Move read_size metadata from the internal namespace to readonly.

void migrate_003_004(manifest& dataset)
{
	auto& metadata = dataset.metadata();
	if (metadata.has_key(meta::ns::internal, "read_size")) {
		uint64_t read_size = metadata.get<uint64_t>(meta::ns::internal, "read_size");
		metadata.set(meta::ns::readonly, "read_size", read_size);
		metadata.unset(meta::ns::internal, "read_size");
		SPLOG("migrate_003_004> (%s, %s, %lu) -> (%s, %s, %lu)",
			meta::ns::internal,
			"read_size",
			read_size,
			meta::ns::readonly,
			"read_size",
			read_size
		);
	}
}

#include "modules/mapred/migrate.h"
#include "modules/io/log.h" 

namespace 
{

typedef std::pair<std::string, std::string> ns_key_t;
typedef std::map<std::string, ns_key_t> table_t;

const table_t conversion {
	{ "kmer_size", { meta::ns::user,     "kmer_size" } },
	{ "kmer_db",   { meta::ns::internal, "kmer_db"   } },
	{ "read_size", { meta::ns::internal, "read_size" } },
	{ "encoding",  { meta::ns::internal, "encoding"  } },
	{ "entries",   { meta::ns::internal, "entries"   } }
};

} // namespace (anonymous)

//
// Move metadata from "tags" into "all_metadata"
//
void migrate_000_001(manifest& dataset)
{
	if (dataset.tags.empty()) {
		SPLOG("migrate_000_001> no (key,value) pairs to be migrated");
		return;
	}

	for (auto tag : dataset.tags) {
		auto it = conversion.find(tag.first);
		if (it != conversion.cend()) {
			if (tag.first == "encoding") {
				// special case for encoding, old tag is a json-serialized string.
				std::string value;
				json_deserialize(value, tag.second);
				dataset.metadata().set(it->second.first, it->second.second, value);
			}
			else {
				dataset.metadata().set(it->second.first, it->second.second, tag.second);
			}
			SPLOG("migrate_000_001> (%s, <value>) -> (%s, %s, <value>)",
				tag.first.c_str(),
				it->second.first.c_str(),
				it->second.second.c_str()
			);
		}
		else {
			SPLOG("migrate_000_001> deleting key: %s", tag.first.c_str());
		}
	}
	dataset.tags.clear();
}

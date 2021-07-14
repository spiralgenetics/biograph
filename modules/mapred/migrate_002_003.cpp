#include "modules/mapred/migrate.h"
#include "modules/io/log.h" 

namespace 
{

typedef std::pair<std::string, std::string> ns_key_t;
typedef std::set<ns_key_t> table_t;

const table_t conversions {
	{ meta::ns::internal, "created" },
	{ meta::ns::internal, "read_size" },
	{ meta::ns::internal, "entries" },
	{ meta::ns::readonly, "kmer_size" },
	{ meta::ns::readonly, "sample_bases" },
	{ meta::ns::readonly, "corrected_read_count" },
	{ meta::ns::readonly, "corrected_base_dist" },
	{ meta::ns::readonly, "failed_correction_count" },
};

} // namespace (anonymous)

//
// Convert JSON encoded strings to integer values
//
void migrate_002_003(manifest& dataset)
{
	auto& metadata = dataset.metadata();
	for (const auto& ns_key : conversions) {
		auto ns = ns_key.first;
		auto key = ns_key.second;
		if (metadata.has_key(ns, key)) {
			std::string json = metadata.get<std::string>(ns, key);
			if (key == "corrected_base_dist") {
				std::vector<uint64_t> value;
				json_deserialize(value, json);
				metadata.set(ns, key, value);
			}
			else {
				size_t value;
				json_deserialize(value, json);
				metadata.set(ns, key, value);
			}
			SPLOG("migrate_002_003> converting (%s, %s, %s)",
				ns.c_str(),
				key.c_str(),
				json.c_str()
			);
		}
	}
}

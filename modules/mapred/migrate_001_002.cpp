#include "modules/mapred/migrate.h"
#include "modules/io/log.h" 

namespace 
{

typedef std::pair<std::string, std::string> ns_key_t;
typedef std::map<ns_key_t, std::string> table_t;

const table_t conversion {
	{ { meta::ns::user, "kmer_size"}, meta::ns::readonly },
	{ { meta::ns::user, "version"},   meta::ns::readonly },
	{ { meta::ns::user, "created"},   meta::ns::readonly }
};

} // namespace (anonymous)

//
// Move the following metadata to the "spiral_readonly" namespace
//   spiral/kmer_size
//   spiral/version
//   spiral/created
//
void migrate_001_002(manifest& dataset)
{
	auto& metadata = dataset.metadata();
	for (const auto& ns_k : conversion) {
		if (metadata.has_key(ns_k.first.first, ns_k.first.second)) {
			std::string value = metadata.get<std::string>(ns_k.first.first, ns_k.first.second);
			metadata.set(ns_k.second, ns_k.first.second, value);
			metadata.unset(ns_k.first.first, ns_k.first.second);
			SPLOG("migrate_001_002> (%s, %s, %s) -> (%s, %s, %s)",
				ns_k.first.first.c_str(),
				ns_k.first.second.c_str(),
				value.c_str(),
				ns_k.second.c_str(),
				ns_k.first.second.c_str(),
				value.c_str()
			);
		}
		else {
			SPLOG("migrate_001_002> (%s, %s) was not found in this dataset ",
				ns_k.first.first.c_str(),
				ns_k.first.second.c_str()
			);
		}
	}
}

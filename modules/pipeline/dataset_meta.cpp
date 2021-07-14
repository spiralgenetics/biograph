#include "modules/pipeline/dataset_meta.h"
#include "modules/pipeline/dataset_path.h"
#include "modules/pipeline/ottoman.h"
#include "modules/web/couchdb.h"
#include "modules/mapred/path.h"
#include "modules/io/config.h"

void dataset_meta::validate() {
	// SPLOG_P(LOG_DEBUG, "dataset_meta::validate> ref_name: %s, the_manifest: %lu, sort_keys: %lu, in_progress: %s",
	// 	ref_name.c_str(),
	// 	the_manifest.get_size(),
	// 	sort_keys.size(),
	// 	in_progress ? "true" : "false"
	// );
}

void dataset_meta::merge_add_element(dataset_meta& total, const std::string& part_url)
{
	SPLOG("dataset_meta::merge_add_element> %s", part_url.c_str());

	couch_server<direntry> db(ottoman_url());
	dataset_path part_path(part_url);
	path::exist_enum e = part_path.exists();
	if (e == path::e_no_exist) {
		throw io_exception("File not found during merge");
	}
	if (e == path::e_directory) {
		SPLOG("dataset_meta::merge_add_element> directory");
		auto out = db.find_match<direntry>("by_parent", part_path.url());
		for (const auto& part : out) {
			merge_add_element(total, part_path.append(part.name).url());
		}
	}
	else {
		dataset_meta meta;
		part_path.load(meta);
		if (!total.type) {
			total.type = meta.type;
		}
		else if (total.type != meta.type) {
			throw io_exception("Data types of merge do not match");
		}

		if (!meta.sort_keys.empty()) {
			throw io_exception("Can't merge sorted files");
		}

		total.the_manifest.add(meta.the_manifest);
	}
}

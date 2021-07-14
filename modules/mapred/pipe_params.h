#pragma once

#include <vector>
#include <string>
#include <map>
#include <limits>

#include "modules/io/transfer_object.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/temp_file.h"

struct pipe_params
{
	TRANSFER_OBJECT {
		VERSION(0);

		FIELD(command, TF_STRICT);
		FIELD(args);
		FIELD(working_dir);
		FIELD(importer_type, TF_STRICT);
		FIELD(exporter_type, TF_STRICT);
		FIELD(ex_im_porter_data);
		FIELD(temp_files);
		FIELD(are_keys_sorted);
	}

	pipe_params() : are_keys_sorted(false) {}

	std::string command;
	std::vector<std::string> args;
	std::string working_dir;
	std::string importer_type;
	std::string exporter_type;
	std::string ex_im_porter_data;
	std::vector<temp_file_spec> temp_files;
	bool are_keys_sorted;
};

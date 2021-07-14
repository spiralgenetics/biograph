#pragma once

#include "modules/io/transfer_object.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/sorter.h"

class kv_sink;

// NOTE: if the 'encoding' field is left empty before calling 'build'
// then output_stream_params will override it and set 'encoding' to 'gzip'
// If you do not want to have any encoding applied to the produced chunks,
// you need to set 'encoding' to 'null'.
//
class output_stream_params
{ 
public:
	TRANSFER_OBJECT
	{ 
		VERSION(0);
		FIELD(unique_str, TF_STRICT); 
		FIELD(goal_size, TF_STRICT); 
		FIELD(num_partitions, TF_STRICT); 
		FIELD(presorted, TF_STRICT); 
		FIELD(allow_split, TF_STRICT); 
		FIELD(clean_break, TF_STRICT);
		FIELD(begin_on, TF_STRICT);
		FIELD(end_before, TF_STRICT);
		FIELD(sort, TF_STRICT); 
		FIELD(reduce, TF_STRICT); 
		FIELD(reduce_param, TF_STRICT); 
		FIELD(split, TF_STRICT); 
		FIELD(encoding);
	}

	std::string unique_str;
	size_t goal_size = 64*1024*1024;
	size_t num_partitions = 1;
	bool presorted = false;  // Used when input will already be in sorted order
	bool allow_split = false;
	bool clean_break = false;
	std::string begin_on; // Start at the first clean record boundry >= here, if empty, use all
	std::string end_before; // Go as long as we are < here, if empty, use all
	std::string sort;
	std::string split;
	std::string reduce; // Implies 'summarize'
	std::string reduce_param;
	std::string encoding;

	std::unique_ptr<kv_sink> build(
		const path& base_path, 
		const std::string& name_prefix, 
		manifest& out
	);
};

#pragma once

#include "modules/mapred/merge_reader.h"
#include "modules/mapred/multi_reader.h"

#include <vector>

class input_stream_params
{ 
public:
	TRANSFER_OBJECT
	{ 
		VERSION(0);
		FIELD(num_records, TF_STRICT); 
		FIELD(inputs, TF_STRICT); 
		FIELD(sort, TF_STRICT); 
		FIELD(clean_break, TF_STRICT);
		FIELD(begin_on, TF_STRICT);
		FIELD(end_before, TF_STRICT); 
		FIELD(split, TF_ALLOW_NULL); 
		FIELD(encoding);
	}

	size_t num_records = 0;	
	std::vector<file_info> inputs; 
	std::string sort;
	std::string split;
	bool clean_break = false;
	std::string begin_on; // Start at this key (or this group if 'clean_break'), ignored if empty
	std::string end_before; // End before this key (or this group if 'clean_break'), ignored if empty
	std::string encoding; // how to decode the underlying file_info data, usually via gzip

	std::unique_ptr<kv_source> build()
	{
		if (sort != "") {
			return std::unique_ptr<kv_source>(
				new merge_reader(sort, 
					inputs.begin(), inputs.end(), 
					begin_on, end_before, 
					clean_break, 
					encoding
				)
			);
		}
		return std::unique_ptr<kv_source>(make_multi_reader(inputs.begin(), inputs.end(), encoding));
	}
};

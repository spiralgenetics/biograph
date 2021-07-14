#pragma once

#include "modules/mapred/manifest.h"
#include "modules/io/progress.h"
#include "modules/io/parallel.h"

// Function that parallelizes a functor/lambda/function ptr over the file_infos
// in a manifest.  The function must have the signature
// void f(const key_type&, const value_type&, size_t file_info_id, size_t record_id)
//
// Where file_info_id is a unique ID between 0 and file_info count - 1 telling
// the function which file_info it has and record_id is the number of the record
// between 0 and the number of KV records in the file_info being processed.
//
// Because the function is called in parallel, be very careful about maintaining
// non-const global state in the function.  If you want mutable state, you probably
// want a vector of state indexed by file_info_id.

template<typename Function, typename KeyType, typename ValueType>
inline Function manifest_parallelize(manifest the_manifest, Function f, progress_handler_t progress = null_progress_handler)
{
	std::vector<file_info> file_infos;
	std::vector<size_t> cumulative_record_ids(1);
	for (auto the_file_info : the_manifest) {
		file_infos.push_back(the_file_info);
		cumulative_record_ids.push_back(cumulative_record_ids.back() + the_file_info.num_records);
	}

    parallel_for(0, file_infos.size(), [&](size_t i) {
		auto file_reader_ptr = file_infos[i].file.read();
		auto file_decoder = make_decoder(the_manifest.get_encoding(), *file_reader_ptr);
		kv_reader file_kv_reader(*file_decoder);
		size_t j = 0;
		KeyType key;
		ValueType value;
		while (file_kv_reader.read_msgpack(key, value)) {
			f(key, value, i, cumulative_record_ids[i] + j++);
		}
      }, progress);
	
	return f;
}

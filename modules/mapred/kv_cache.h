#ifndef __kv_cache_h__
#define __kv_cache_h__

#include "modules/mapred/kv_hold.h"

// kv_cache is a rw cache for key-value pairs.
// It is not meant to be used directly, but rather via a base_chunker instance.
// The base chunker will call kv_cache::oversized to check whether a chunk of data
// needs to be written to a file.
//
// Example:
// #include "mapred/manifest.h"
// #include "mapred/base_chunker.h"
// #include "mapred/kv_cache.h"
//
// manifest out;
// path manifest_path("out/test_out/");
// simple_name_generator fng(manifest_path.append("input"));
// size_t records_per_chunk = 1000*1000;
//
// base_chunker<kv_cache> out_chunker("lexical", fng, "chunk", records_per_chunk, 0, out);
//
// // generate a 1000 "chunk" files, each with 1 million values in them
// for(size_t i=0; i< 1000*1000*1000; i++)
// {
//   key = printstring("%15ld",i);
//   out_chunker.write(key, key); //value = key in this example.
// }

class kv_cache : public kv_hold
{
public:
	kv_cache(const std::string& sort):kv_hold(sort) {}
	bool oversized(size_t records_per_chunk) const { return get_num_records() >= records_per_chunk; } 
};

#endif

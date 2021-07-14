#include "modules/test/test_utils.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/reduce_task.h"
#include "modules/mapred/sort_task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/map_task.h"
#include "modules/mapred/kv_cache.h"
#include "modules/mapred/base_chunker.h"
#include "modules/io/config.h"
#include "modules/io/log.h"
#include "gtest/gtest.h"

#include <map>

// generate a manifest full of chunks of files. Each file has ordered data, and overlaps with the previous and next chunk, by one key.
// For example:
// chunk_a  0, 1, 2, ... 98,     100 <EOF>
// chunk_b                   99,     101, 102, ..., 198,      200 <EOF>
// chunk_c                                               199,     201, 202, ... <EOF>
void gen_kv(const path& manifest_path, const size_t numKV, const size_t records_per_chunk,  manifest& out)
{
	base_chunker<kv_cache> out_chunker("lexical", manifest_path.append("input"), "chunk", records_per_chunk, 0, out, codec::gzip);
	std::string key,key_tmp;

	key.resize(20);
	key_tmp.resize(20);
	size_t remainder, num_1(records_per_chunk -1);
	size_t i = 0;
	while( i < numKV)
	{
		remainder = i % records_per_chunk;
		if( remainder == num_1)
		{
			key_tmp = printstring("%12ld",i);
			i++;
			key = printstring("%12ld",i);
			out_chunker.write(key, key); //value = key
			out_chunker.write(key_tmp, key_tmp); //value = key
		}
		else
		{
			key = printstring("%12ld",i);
			out_chunker.write(key, key); //value = key
		}
		i++;
	}
	out_chunker.close();
}

TEST(merge_sort, empty)
{
	//make sure split-sorting an empty manifest returns an empty vector
	std::vector<input_stream_params> inputs;
	manifest empty("lexical");
	empty.split_sort(empty, inputs, 30);
	ASSERT_TRUE(inputs.empty());
}

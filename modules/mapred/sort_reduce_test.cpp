#include "base/base.h"
#include "modules/test/test_utils.h"
#include "gtest/gtest.h"
#include "modules/test/test_utils.h"
#include "modules/mapred/sort_task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/map_task.h"
#include "modules/mapred/base_chunker.h"
#include "modules/mapred/kv_hold.h"
#include "modules/mapred/reducer.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"
#include "modules/bio_base/struct_var.h"

class test_reducer : public typed_reducer<test_reducer, struct_var_key, int, int, int>
{
public:
        test_reducer(const std::string& params) {}

        void typed_start(const struct_var_key& key) { m_key = key.variation_id, m_value = 0; }
        void typed_add_value(const struct_var_key& key, const int& value) { m_value += value; }
        void typed_end() { output(m_key * 5123 % 1049, m_value); } 

private:
	int m_key;
	int m_value;
};

REGISTER_1(reducer, test, const std::string&);

TEST(sort_reduce_test, gzip)
{
	manifest orig_manifest;
	path test_path(make_path("sort_test"));
	std::map<struct_var_key, int> verify; 

        base_chunker<kv_hold> out_chunker("", test_path.append("input"), "chunk", 10*1024*1024, 0, orig_manifest, codec::gzip);
	int next_read = 0;
        for (size_t i = 0; i < 10000; i++) {
		struct_var_key k;
		k.variation_id = random() % 1000;
		k.read_id = next_read++;
		int v = random() % 1000;
		verify[k] = v;
		out_chunker.write_msgpack(k, v);
	}
        out_chunker.close();

	size_t chunk_size = 5000; // Room for about 100 values

	task_mgr_local tm;
	std::unique_ptr<map_task> mt = make_unique<map_task>();
	mt->input = orig_manifest;
	mt->map = "identity";
	mt->output_goal_size = chunk_size;
	mt->sort = "struct_var";

	SPLOG("--- Running do_map 'identity'");
	manifest map_manifest;
	tm.run_task(map_manifest, test_path.append("do_map"), std::move(mt));

	SPLOG("--- Done with do_map 'identity'");

	std::unique_ptr<sort_task> st = make_unique<sort_task>();
	st->input = map_manifest;
	st->goal_size = chunk_size;
	st->max_files = 8;

	SPLOG("--- Running do_sort sort_task");
	manifest sort_manifest;
	tm.run_task(sort_manifest, test_path.append("do_sort"), std::move(st));
	SPLOG("--- Done with  do_sort sort_task");

	SPLOG("Sort manifest records = %d", (int) sort_manifest.get_num_records());
	SPLOG("Manifest as string: %s", json_serialize(sort_manifest).c_str());
	auto it = verify.begin();
	manifest_reader kv_read(sort_manifest);

	struct_var_key key;	
	int value;
	unsigned int cur_var_id = it->first.variation_id;
	int total = 0;
	std::map<int, int> reduced;
	while (kv_read.read_msgpack(key, value)) {
		if (cur_var_id != key.variation_id) {
			reduced[cur_var_id * 5123 % 1049] = total;
			//reduced[cur_var_id] = total;
			cur_var_id = key.variation_id;
			total = 0;
		}
		total += it->second;
		//SPLOG("key = %s, value = %s", json_serialize(key).c_str(), json_serialize(value).c_str());
		ASSERT_EQ(key, it->first);
		ASSERT_EQ(value, it->second);
		++it;
	}
	CHECK(it == verify.end());
	reduced[cur_var_id * 5123 % 1049] = total;

        std::unique_ptr<sorted_reduce_task> srt = make_unique<sorted_reduce_task>();
        srt->input = sort_manifest;
	srt->reduce = "test";
	srt->out_sort = "lexical";
	srt->prereduce_goal_size = chunk_size;

	SPLOG("--- Running sorted reduce task"); 
	manifest reduced_manifest;
        tm.run_task(reduced_manifest, test_path.append("do_reduce"), std::move(srt));
	SPLOG("--- Done with sorted reduce task");

	st = make_unique<sort_task>();
	st->input = reduced_manifest;

	SPLOG("--- Running final sort task"); 
	manifest reduced_2_manifest;
        tm.run_task(reduced_2_manifest, test_path.append("do_xsort"), std::move(st));
	SPLOG("--- Done with sorted reduce task");

	int ikey;
	auto it2 = reduced.begin();
	manifest_reader kv_read2(reduced_2_manifest);
	while (kv_read2.read_msgpack(ikey, value)) {
		//SPLOG("key = %d, %d value = %d, %d", ikey, it2->first, value, it2->second);
		ASSERT_EQ(ikey, it2->first);
		ASSERT_EQ(value, it2->second);
		++it2;
	}
	CHECK(it == verify.end());
}


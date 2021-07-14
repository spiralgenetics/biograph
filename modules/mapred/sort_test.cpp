#include "base/base.h"
#include "modules/test/test_utils.h"
#include "gtest/gtest.h"
#include "modules/test/test_utils.h"
#include "modules/mapred/sort_task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/map_task.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"

void main_sort_test(const std::string& encoding)
{
	manifest orig_manifest;
	path test_path(make_path("sort_test"));
	std::map<std::string, std::string> verify;
	gen_random_kv(test_path, 10000, 10*1024*1024, 20, verify, orig_manifest, encoding);
	size_t chunk_size = 5000; // Room for about 100 values

	task_mgr_local tm;
	std::unique_ptr<map_task> mt = make_unique<map_task>();
	mt->input = orig_manifest;
	mt->map = "identity";
	mt->output_goal_size = chunk_size;
	mt->sort = "lexical";

	SPLOG("--- Running do_map 'identity'");
	manifest map_manifest;
	tm.run_task(map_manifest, test_path.append("do_map"), std::move(mt));

	SPLOG("--- Done with do_map 'identity'");

	ASSERT_EQ(map_manifest.metadata().get<std::string>(meta::ns::internal, "encoding"), encoding);

	std::unique_ptr<sort_task> st = make_unique<sort_task>();
	st->input = map_manifest;
	st->goal_size = chunk_size;
	st->max_files = 8;

	SPLOG("--- Running do_sort sort_task");
	manifest sort_manifest;
	tm.run_task(sort_manifest, test_path.append("do_sort"), std::move(st));
	SPLOG("--- Done with  do_sort sort_task");

	ASSERT_EQ(sort_manifest.metadata().get<std::string>(meta::ns::internal, "encoding"), encoding);

	SPLOG("Sort manifest records = %d", (int) sort_manifest.get_num_records());
	SPLOG("Manifest as string: %s", json_serialize(sort_manifest).c_str());
	std::map<std::string, std::string>::const_iterator it = verify.begin();
	manifest_reader kv_read(sort_manifest);
	std::string key,value;
	while (kv_read.read(key, value)) {
		//SPLOG("key = %s, value = %s", key.c_str(), value.c_str());
		ASSERT_EQ(key, it->first);
		ASSERT_EQ(value, it->second);
		++it;
	}
	CHECK(it == verify.end());
}

TEST(sort_test, gzip)
{
	main_sort_test(codec::gzip);
}

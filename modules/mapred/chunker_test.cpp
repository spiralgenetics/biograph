#include "modules/test/test_utils.h"
#include "modules/mapred/base_chunker.h"
#include "modules/mapred/kv_hold.h"
#include <gtest/gtest.h>

void generate_and_compare(const std::string& encoding)
{
	kv_hold compare("");
	path test_path(make_path("chunker"));
	manifest chunk_manifest;
	base_chunker<kv_hold> out("", test_path, "chunk", 100, 0, chunk_manifest, encoding);

	std::string key, value;
	key.resize(20);
	value.resize(20);
	for(size_t i = 0; i < 20; i++)
	{
		for(size_t i = 0; i < 20; i++)
		{
			key[i] = 'a' + random()%26;
			value[i] = 'a' + random()%26;
		}
		out.write(key, value);
		compare.write(key, value);
	}
	out.close();

	manifest_reader readback(chunk_manifest);
	while(readback.read(key, value))
	{
		std::string key2, value2;
		bool r = compare.read(key2, value2);
		ASSERT_TRUE(r);
		ASSERT_EQ(key, key2);
		ASSERT_EQ(value, value2);
	}
	ASSERT_FALSE(compare.read(key, value));
}

TEST(chunker, null)
{
	generate_and_compare(codec::null);
}

TEST(chunker, gzip)
{
	generate_and_compare(codec::gzip);
}

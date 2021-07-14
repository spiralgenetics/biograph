#include "modules/test/test_utils.h"
#include "modules/io/keyvalue.h"
#include "modules/mapred/kv_hold.h"
#include "modules/io/inverse_kvread.h"
#include "modules/mapred/path.h"
#include <gtest/gtest.h>

void make_random(std::string& s)
{
	size_t size = (2 << (random() % 15)) + random() % 200;
	s.clear();
	for(size_t i = 0; i < size; i++)
		s.push_back(rand() % 26 + 'a');
}

void make_random(std::string& key, std::string& value)
{
	make_random(key);
	make_random(value);
}

TEST(keyvalue, normal)
{
	kv_hold cmp("");
	path file(make_path("keyvalue_normal"));
	std::unique_ptr<writable> writer = file.write();
	kv_writer wx(*writer);
	std::string key, value;
	for(size_t i = 0; i < 1000; i++)
	{
		make_random(key, value);
		cmp.write(key, value);
		wx.write(key, value);
	}
	wx.close();
	writer->close();
	std::unique_ptr<readable> reader = file.read();
	kv_reader rx(*reader);
	while(rx.read(key, value))
	{
		std::string key2, value2;
		bool r = cmp.read(key2, value2);
		ASSERT_TRUE(r);
		ASSERT_EQ(key, key2);
		ASSERT_EQ(value, value2);
	}
	ASSERT_FALSE(cmp.read(key, value));
}

TEST(keyvalue, inverse)
{
	kv_hold cmp("");
	path file(make_path("keyvalue_inverse"));
	std::string key, value;
	for(size_t i = 0; i < 1000; i++)
	{
		make_random(key, value);
		cmp.write(key, value);
	}
	inverse_kvread kvr(cmp);
	std::unique_ptr<writable> writer = file.write();
	io_copy(kvr, *writer);
	writer->close();
	cmp.reset();
	std::unique_ptr<readable> reader = file.read();
	kv_reader rx(*reader);
	while(rx.read(key, value))
	{
		std::string key2, value2;
		bool r = cmp.read(key2, value2);
		ASSERT_TRUE(r);
		ASSERT_EQ(key, key2);
		ASSERT_EQ(value, value2);
	}
	ASSERT_FALSE(cmp.read(key, value));
}



#include "modules/mapred/metadata.h" 
#include "gtest/gtest.h"
#include "modules/io/log.h"

meta::merge::init registration(meta::merge::register_fn("metadata_unittest_e", [](
		const meta::merge::params& params) 
	{ 
		return params.value2; 
	} 
));

TEST(metadata, basic)
{
	meta::data m;

	int a_value = 13;
	m.set("foo", "bar", a_value);
	ASSERT_TRUE(m.has_key("foo", "bar"));
	ASSERT_EQ(m.get<int>("foo", "bar", 0), a_value);

	ASSERT_EQ(m.get<int>("foo", "bar"), a_value);
	ASSERT_EQ(m.get("foo", "yo", 47), 47);

	ASSERT_THROW(m.get<int>("", "bar"), io_exception);
	ASSERT_THROW(m.get<int>("foo", ""), io_exception);

	// this key does not exist
	ASSERT_FALSE(m.has_key("foo", "bad_key"));
	ASSERT_EQ(m.get<int>("foo", "bad_key", 11),	11);
	ASSERT_THROW(m.get<int>("foo", "bad_key"), io_exception);

	// this namespace does not exist
	ASSERT_EQ(m.get<int>("bad_ns", "bar", 253),	253);
	ASSERT_THROW(m.get<int>("bad_ns", "bar"), io_exception);
}

TEST(metadata, merge)
{
	meta::data m1, m2;

	m1.set("a", "b", std::string("c"));
	m2.set("a", "b", std::string("c"));
	m1.set("x", "b", std::string("c"));
	m2.set("y", "b", std::string("c"));
	m1.set("n", "k", std::string("c"));
	m2.set("n", "K", std::string("C"));

	ASSERT_NO_THROW(m1.merge(m2));
	ASSERT_EQ(m1.get<std::string>("a", "b"), std::string("c"));
	ASSERT_EQ(m1.get<std::string>("x", "b"), std::string("c"));
	ASSERT_EQ(m1.get<std::string>("y", "b"), std::string("c"));
	ASSERT_EQ(m1.get<std::string>("n", "k"), std::string("c"));
	ASSERT_EQ(m1.get<std::string>("n", "K"), std::string("C"));

	meta::data m3,m4;
	//collision
	m3.set("a", "d", std::string("e"));
	m4.set("a", "d", std::string("f"));
	ASSERT_THROW(m3.merge(m4), io_exception);

	//try a custom collision handler
	meta::data m5,m6;
	m5.set("a", "metadata_unittest_e", std::string(json_serialize(std::string("e"))));
	m6.set("a", "metadata_unittest_e", std::string(json_serialize(std::string("f"))));
	ASSERT_NO_THROW(m5.merge(m6));
	ASSERT_EQ(m5.get<std::string>("a", "metadata_unittest_e"), json_serialize(std::string("f")));
}

TEST( metadata, unset )
{
	meta::data m;

	m.set("foo", "bar", std::string("joe"));
	ASSERT_EQ(m.get<std::string>("foo", "bar"), std::string("joe"));
	ASSERT_NO_THROW(m.unset("foo", "bar"));
	ASSERT_THROW(m.get<std::string>("foo", "bar"), io_exception);
	ASSERT_NO_THROW(m.unset("foo", "bad_key"));
}


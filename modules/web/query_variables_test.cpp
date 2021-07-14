#include "gtest/gtest.h"
#include "modules/web/url_query.h" 

TEST(qv, a)
{
	std::string qs1("name1=value1&name2=value2&name3=value3");
	query_variables vars_1(qs1),
			vars_2("paired=false");

	EXPECT_NO_THROW( query_variables vars_3("somecrap"));

	query_variables e_vars_1 = std::map<std::string,std::string> {
		{"name1", "value1"},
		{"name2", "value2"},
		{"name3", "value3"}
	};

	query_variables e_vars_2 = std::map<std::string,std::string> {
		{"paired", "false"}
	};

	ASSERT_TRUE( vars_1 == e_vars_1 );
	ASSERT_TRUE( vars_2 == e_vars_2 );
	ASSERT_EQ(0, qs1.compare(vars_1));

	query_variables empty_qv(""), empty_too("&");
	ASSERT_TRUE(empty_qv.empty());
	ASSERT_TRUE(empty_too.empty());

	query_variables single_pair("name=value");
	ASSERT_EQ( (size_t)1, single_pair.size());
	ASSERT_EQ(0, std::string("name=value").compare(single_pair));

	ASSERT_EQ(0, std::string("name=").compare(query_variables("name=")));
}

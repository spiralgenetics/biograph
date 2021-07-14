#include "gtest/gtest.h"
#include "modules/web/urlencode.h"

TEST(urlencode, url)
{
	std::string url("a/ c/{e/f}");
	ASSERT_EQ( urlencode(url).compare("a%2F%20c%2F%7Be%2Ff%7D"), 0);
	ASSERT_EQ( urldecode(urlencode(url)).compare(url), 0);
}

TEST(urlencode, component)
{
	std::string url("a/ c/{e/f}");
	ASSERT_EQ( urlencode_component(url).compare("a/%20c/%7Be/f%7D"), 0);
	ASSERT_EQ( urldecode(urlencode_component(url)).compare(url), 0);

	std::string no_forward_slash("a c{>");
	ASSERT_EQ( urlencode_component( no_forward_slash ).compare("a%20c%7B%3E"), 0);

	std::string contiguous_forward_slashes(" // ");
	ASSERT_EQ( urlencode_component( contiguous_forward_slashes ).compare("%20/%20"), 0); // notice the effect of boost::token_compress_on

	std::string leading_and_trailing("/ /");
	ASSERT_EQ( urlencode_component( leading_and_trailing ).compare("/%20/"), 0);
}

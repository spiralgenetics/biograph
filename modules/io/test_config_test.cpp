
#include "gtest/gtest.h"
#include "modules/io/config.h"
#include "modules/io/io.h"
#include <stdio.h>
#include <boost/filesystem.hpp>

TEST(Config, ValidateDefaultConfig)
{
	boost::filesystem::path path("etc/products/unittest.json");
	ASSERT_NO_THROW(Config::load(path.string()));
	EXPECT_FALSE(CONF_T(std::string, url_base).empty());
	EXPECT_FALSE(CONF_T(std::string, path_bulkdata).empty());
	EXPECT_FALSE(CONF_T(std::string, reference_path).empty());
	EXPECT_FALSE(CONF_T(std::string, path_user_base).empty());
	EXPECT_FALSE(CONF_T(std::string, path_reference_base).empty());
}

TEST(Config, InvalidConfig)
{
	//file that does not exist
	EXPECT_THROW({ Config::load("/wheredoyouthinkyouare?"); }, io_exception);
	//Throw on valid file with invalid JSON
	EXPECT_THROW({ Config::load("config/invalid-config.json"); }, io_exception);
}

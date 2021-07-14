#include "gtest/gtest.h"
#include "modules/io/version.h"
#include "modules/io/log.h"

void test_version(
	const std::string& version,
	unsigned major,
	unsigned minor,
	unsigned patch,
	const std::string& pre,
	const std::string& build)
{
	SPLOG("Testing version: %s", version.c_str());
	product_version vp(version);
	ASSERT_EQ(major,     vp.major);
	ASSERT_EQ(minor,     vp.minor);
	ASSERT_EQ(patch,     vp.patch);
	ASSERT_EQ(pre,       vp.pre);
	ASSERT_EQ(build,     vp.build);
	ASSERT_EQ(version,   vp.make_string());
}

TEST(version, helper)
{
	test_version("1.2.3", 1, 2, 3, "", "");
	test_version("1.2.3-pre.1", 1, 2, 3, "pre.1", "");
	test_version("1.2.3+This-is-a.build", 1, 2, 3, "", "This-is-a.build");
	test_version("1.2.3-pre.1+This-is-a.build", 1, 2, 3, "pre.1", "This-is-a.build");
	ASSERT_THROW(product_version("1.2"), io_exception);
	ASSERT_THROW(product_version("1.2.3-pre-pre"), io_exception);
	ASSERT_THROW(product_version("1.2.3+X!"), io_exception);
}


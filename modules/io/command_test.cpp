#include "modules/io/command.h"
#include "modules/io/mem_io.h"
#include "modules/io/log.h"
#include <gtest/gtest.h>
#include <boost/algorithm/string.hpp>

TEST(command, basic)
{
  mem_io input("", track_alloc("command_test"));
	input.print("Hello\nWorld\n");
	auto output = exec::communicate(input, "grep", {"World"});
	ASSERT_EQ("World\n", output);
}

TEST(command, exe_dir) 
{
	std::string exe_dir = exec::get_exe_dir(::getpid());
	SPLOG("This process is running from directory %s", exe_dir.c_str());

	EXPECT_TRUE(boost::algorithm::ends_with(exe_dir, "modules/io"));
}

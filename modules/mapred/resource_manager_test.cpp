#include "modules/test/test_utils.h"
#include "modules/mapred/resource_manager.h"
#include "modules/io/json_transfer.h"
#include <gtest/gtest.h>

TEST(resource_manager, direct)
{
	resource_manager rm(true);
	mmap_buffer mmb;
	rm.create_resource(mmb, 1024);

	std::string rpath = mmb.path();
	SPLOG("Path = %s", rpath.c_str());
	strcpy(mmb.buffer(), "Hello World");

	path root(make_path("rm"));
	manifest out;
	SPLOG("Writing resource");
	rm.write_resource(out, mmb, root, "prefix"); 
	SPLOG("Manifest = %s", json_serialize(out).c_str());
	mmb.close();

	rm.read_resource(mmb, out);
	ASSERT_EQ(strcmp(mmb.buffer(), "Hello World"), 0);
	mmb.close();
}

TEST(resource_manager, indirect)
{
	resource_manager rm(false);
	mmap_buffer mmb;
	rm.create_resource(mmb, 1024);

	std::string rpath = mmb.path();
	SPLOG("Path = %s", rpath.c_str());
	strcpy(mmb.buffer(), "Hello World");

	path root(make_path("rm"));
	manifest out;
	SPLOG("Writing resource");
	rm.write_resource(out, mmb, root, "prefix"); 
	SPLOG("Manifest = %s", json_serialize(out).c_str());
	mmb.close();

	rm.read_resource(mmb, out);
	ASSERT_EQ(strcmp(mmb.buffer(), "Hello World"), 0);
	mmb.close();

	unlink(rpath.c_str());
	rm.read_resource(mmb, out);
	ASSERT_EQ(strcmp(mmb.buffer(), "Hello World"), 0);
}

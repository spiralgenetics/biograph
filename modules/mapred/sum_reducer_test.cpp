
#include "gtest/gtest.h"
#include "modules/mapred/sum_reducer.h"
#include "modules/io/make_unique.h"

TEST(sum_reducer, summarize)
{
	uint64_t x = 5;
	uint64_t y = 7;

	std::string x1 = msgpack_serialize(x);
	std::string y1 = msgpack_serialize(y);

	std::unique_ptr<reducer> r = make_unique<sum_reducer>("");
	r->summarize(x1, y1);
	uint64_t o = 0;
	msgpack_deserialize(o, x1);
	
	ASSERT_EQ(o, x+y);
}

#include "modules/io/mem_io.h"
#include "modules/io/log.h"
#include "modules/io/readable_PRNG.h"
#include "gtest/gtest.h"

TEST(mem_io, basic)
{
  mem_io mio("", track_alloc("mem_io_test"));
	mio.write("Hello", 5);
	ASSERT_EQ(mio.size(), (unsigned) 5);
	char buf[20];
	size_t r = mio.read(buf, 20);
	ASSERT_EQ(r, (unsigned) 5);
	ASSERT_EQ(memcmp(buf, "Hello", 5), 0);
	mio.reset();
	size_t i;
	for (i = 0; i < 20; i++) {
		r = mio.read(buf + i, 1);
		if (r == 0) {
			break;
		}
	}
	ASSERT_EQ(i, (unsigned) 5);
	ASSERT_EQ(memcmp(buf, "Hello", 5), 0);
}

TEST(mem_io, random)
{
	readable_PRNG prng(16*1024*1024, 4, time(0));
	mem_io buf("", track_alloc("mem_io_test"));

	SPLOG("write random data to read/write buffer");
	io_copy(prng, buf);

	ASSERT_EQ(buf.size(), prng.size());

	SPLOG("read read/write buffer and compare with initial random data");
	prng.reset();
	size_t first_diff_pos;
	ASSERT_TRUE(io_match(buf, prng, first_diff_pos));

	SPLOG("reset read/write buffer and compare again");
	prng.reset();
	buf.reset();
	ASSERT_TRUE(io_match(buf, prng, first_diff_pos));
}

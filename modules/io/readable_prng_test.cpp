#include "modules/io/readable_PRNG.h" 
#include "modules/io/mem_io.h" 
#include "modules/io/zip.h"
#include "modules/io/log.h" 
#include "gtest/gtest.h"

TEST(prng, a)
{
	const size_t size = 64*1024;

	for (size_t randomness = 0; randomness < 9; randomness++) {
		SPLOG("randomness = %zu", randomness);

		readable_PRNG prng(size, randomness, time(0));

		mem_io buf("", track_alloc("readable_prng_test"));
		zip_writer compressor(buf, no_update, 9); // 9 = best zlib compression

		io_copy(prng, compressor);
		compressor.close();

		SPLOG("%3.1f of original size", 100 * buf.size() / double(size));

		if (randomness < 8) {
			ASSERT_LT(buf.size(), size);
		}
	}
}

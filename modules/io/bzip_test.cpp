
#include <gtest/gtest.h>
#include "modules/io/file_io.h"
#include "modules/io/bzip.h"

TEST(bzip, basic)
{
	file_reader fr("golden/test.qseq.bz2");
	file_reader golden_file("golden/test.qseq");
	bzip_reader br(fr);
	size_t first_diff_pos;
	ASSERT_TRUE( io_match( br, golden_file, first_diff_pos ) );
}


#include <stdlib.h>
#include "modules/io/file_io.h"
#include "modules/io/mem_io.h"
#include "modules/io/keyvalue.h"
#include "modules/bio_format/qseq.h"
#include "modules/bio_format/read_qual.h"
#include "gtest/gtest.h"

/*
TEST(QSeq, Import)
{
	file_reader in("/src/golden/test.qseq");

	mem_io buf1;
	mem_io buf2;
	kv_writer buf1w(buf1);
	kv_reader buf1r(buf1);
	kv_writer buf2w(buf2);
	kv_reader buf2r(buf2);

	import_qseq import(in);
	import.import(buf1w);

        job step;
        step.mapper = "read_qual";
        step.reducer = "sum";
        step.sorter = "lexical";

	memory_mapred(step, buf1r, buf2w);
	buf1.buffer.clear();
	file_writer out(stdout);
	export_read_qual export_qual;
	export_qual.do_export(buf2r, out);

	// TODO: Check results!
}
*/

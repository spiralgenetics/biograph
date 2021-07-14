#include "modules/test/test_utils.h"
#include "modules/io/file_io.h"
#include "modules/io/keyvalue.h"
#include "modules/bio_format/fasta.h"
#include "gtest/gtest.h"
#include <stdlib.h>

TEST(FastA, RoundtripFastA)
{
	try {
		file_reader imp_in("golden/sequences.fasta");
		file_writer imp_out(make_path("sequences.kvp"));
		kv_writer kout(imp_out);	
		fasta_importer import(imp_in);
		import.import(kout, discard_simple_metadata());
		imp_in.close();
		imp_out.close();

		file_reader exp_in(make_path("sequences.kvp"));
		auto out_path = make_path("sequences.fasta");
		file_writer exp_out(out_path);
		kv_reader kin(exp_in);	
		fasta_exporter a_fasta_exporter(exp_out);
		a_fasta_exporter.export_from(kin);
		exp_in.close();
		exp_out.close();

		ASSERT_TRUE(diff(out_path, "golden/sequences.fasta"));
	} 
	catch (const io_exception& io) {
		printf("Error: %s\n", io.message().c_str());
		ASSERT_TRUE(false);
	}
}

#include "modules/test/fastq_test_utils.h"
#include "modules/io/file_io.h"
#include "modules/io/keyvalue.h"
#include "modules/io/zip.h"
#include "modules/bio_format/fastq.h"

void make_fastq_kv(const std::string& fastq_file, const std::string& kv_file)
{
	file_reader in(fastq_file);
	file_writer out(kv_file);
	kv_writer kout(out);	
	fastq_importer import(in);
	import.import(kout, discard_simple_metadata());
	in.close();
	out.close();
}


void make_zipped_fastq_kv(const std::string& fastq_file, const std::string& kv_file)
{
	file_reader in(fastq_file);
	zip_reader unzipper(in);
	file_writer out(kv_file);
	kv_writer kout(out);	
	fastq_importer import(unzipper);
	import.import(kout, discard_simple_metadata());
	in.close();
	out.close();
}

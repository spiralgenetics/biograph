#include "modules/io/zip.h"
#include "modules/io/file_io.h"
#include "modules/io/mem_io.h"
#include "modules/io/log.h"
#include "modules/io/readable_PRNG.h" 
#include "gtest/gtest.h"

#include <cstdio> 

void print_update(size_t total, size_t in, size_t out)
{
	//SPLOG("read: %lu, written: %lu", in, out);
	
	double fraction_completed = in/(double)total;

	printf("\33[2K\r%3.0f%%", 100*fraction_completed);
	fflush(stdout);
}

TEST(zip, write)
{
	file_reader test_qseq("golden/test.qseq");

	mem_io buf("", track_alloc("zip_test"));
	size_t input_size = fsize("golden/test.qseq");
	const size_t write_modulo = input_size/1000;
	progress_t update = [&input_size, &write_modulo](size_t in, size_t out) {
		print_update(input_size, in, out);
		return write_modulo;
	};
	zip_writer zipper(buf, update);

	io_copy( test_qseq, zipper );

	// we need the following close: (because we'll re-open the file shortly)
	test_qseq.close();
	zipper.close();
	printf("\n");

	//SPLOG("compressed buffer size: %lu", buf.size());

	zip_reader uncompressed(buf);
	file_reader test_qseq2("golden/test.qseq");

	size_t first_diff_pos;
	ASSERT_TRUE( io_match(uncompressed, test_qseq2, first_diff_pos) );
}

TEST(zip, big)
{
	// STAGE 1: Compressing fake data
	size_t fake_data_size = 240; // MB
	SPLOG("unzip(zip(%luMB))", fake_data_size);

	unsigned int seed = time(0);
	size_t compressible_data_size = fake_data_size*1024*1024;
	readable_PRNG compressible_data(compressible_data_size, 5, seed);
	
	const size_t compressible_modulo = compressible_data_size / 100;
	progress_t update = [&compressible_data_size, &compressible_modulo](size_t in, size_t out) {
		print_update(compressible_data_size, in, out);
		return compressible_modulo;
	};

	mem_io compressed_data("", track_alloc("zip_test"));
	zip_writer zipper(compressed_data, update);

	SPLOG("compressing...");
	io_copy(compressible_data, zipper);

	zipper.close();
	printf("\n");

	compressible_data.reset();

	// STAGE 2: Decompressing compressed fake data
	
	mem_io uncompressed_data("", track_alloc("zip_test"));
	size_t compressed_data_size = compressed_data.size();
	SPLOG("%lu -> %lu", compressible_data_size, compressed_data_size);
	
	const size_t compressed_modulo = compressed_data_size / 100;
	
	progress_t update2 = [&compressed_data_size, &compressed_modulo](size_t in, size_t out) {
		print_update(compressed_data_size, in, out);
		return compressed_modulo;
	};

	zip_reader unzip(compressed_data, update2);

	SPLOG("decompressing...");
	io_copy(unzip, uncompressed_data);
	printf("\n");

	ASSERT_EQ( uncompressed_data.size(), compressible_data.size() );

	// STAGE 3: Verifying decompressed compressed fake data

	SPLOG("verifying uncompressed data...");
	size_t first_diff_pos;

	io_match_update_t match_update = [&compressible_data_size, &compressible_modulo](size_t in) {
		print_update(compressible_data_size, in, 0);
		return compressible_modulo;
	};

	ASSERT_TRUE( io_match(compressible_data, uncompressed_data, first_diff_pos, match_update) );
	printf("\n");
	fflush(stdout);
}

TEST(zip, read)
{
	std::string gz_path("golden/test.qseq.gz");
	std::string qseq_path("golden/test.qseq");

	file_reader test_qseq_gz(gz_path);
	file_reader test_qseq(qseq_path);

	size_t input_size = fsize(gz_path);
	const size_t update_chunk_size = input_size/100;

	progress_t update = [&input_size, &update_chunk_size](size_t in, size_t out) {
		print_update(input_size, in, out);
		return update_chunk_size;
	};

	zip_reader zipper(test_qseq_gz, update);
	mem_io buf("", track_alloc("zip_test"));

	io_copy(zipper, buf);

	printf("\n");

	ASSERT_EQ( (size_t)fsize(qseq_path), buf.size() );
	size_t first_diff_pos;
	ASSERT_TRUE( io_match(test_qseq, buf, first_diff_pos));	
}


TEST(zip, zero_padding)
{
	file_reader test_qseq("golden/test.qseq");

	mem_io buf("", track_alloc("zip_test"));
	size_t input_size = fsize("golden/test.qseq");
	const size_t write_modulo = input_size/1000;
	progress_t update = [&input_size, &write_modulo](size_t in, size_t out) {
		print_update(input_size, in, out);
		return write_modulo;
	};
	zip_writer zipper(buf, update);

	io_copy(test_qseq, zipper);

	// we need the following close: (because we'll re-open the file shortly)
	test_qseq.close();
	zipper.close();
	printf("\n");

	// add zero byte paddding to the end of the compressed buffer
	std::array<char, 100> padding;
	padding.fill(0);
	buf.write(padding.data(), padding.size());

	//SPLOG("compressed buffer size: %lu", buf.size());

	zip_reader uncompressed(buf);
	file_reader test_qseq2("golden/test.qseq");

	size_t first_diff_pos;
	ASSERT_TRUE(io_match(uncompressed, test_qseq2, first_diff_pos));
}


TEST(zip, blocks)
{
	std::string zipped_path("golden/e_coli_10000snp.fq.multiblock.gz");
	std::string unzipped_path("golden/e_coli_10000snp.fq");

	file_reader zipped_file_reader(zipped_path);
	file_reader unzipped_file_reader(unzipped_path);

	size_t input_size = fsize(zipped_path);
	const size_t update_chunk_size = input_size/100;

	progress_t update = [&input_size, &update_chunk_size](size_t in, size_t out) {
		print_update(input_size, in, out);
		return update_chunk_size;
	};

	zip_reader zipper(zipped_file_reader, update);
	mem_io buf("", track_alloc("zip_test"));

	io_copy(zipper, buf);

	printf("\n");

	ASSERT_EQ( (size_t)fsize(unzipped_path), buf.size() );
	size_t first_diff_pos;
	ASSERT_TRUE( io_match(unzipped_file_reader, buf, first_diff_pos));
}

#include "modules/test/test_utils.h"
#include "modules/mapred/output_stream.h"
#include "modules/io/encoding.h"
#include "modules/io/mem_io.h"
#include "modules/io/readable_PRNG.h"
#include "gtest/gtest.h"

#include <memory>

class encoding : public ::testing::Test
{
public:
	encoding()
		: prng(8*1024*1024, 4, time(0))
	{}

protected:
	void SetUp() override { }
	void TearDown() override { }

	readable_PRNG prng;
};

TEST_F(encoding, basic)
{
  mem_io buffer("", track_alloc("encoding_test"));

	ASSERT_THROW(make_encoder("not so random garbage", buffer), unknown_codec);
	ASSERT_THROW(make_decoder("les sanglots longs des violons de l'automne", buffer), unknown_codec);

	std::unique_ptr<writable> compressed_buffer(make_encoder(codec::gzip, buffer));

	io_copy(prng, *compressed_buffer);

	compressed_buffer->close();

	ASSERT_LT(buffer.size(), prng.size());

	std::unique_ptr<readable> uncompressed_buffer(make_decoder(codec::gzip, buffer));

	prng.reset();
	size_t first_diff_pos;
	ASSERT_TRUE(io_match( *uncompressed_buffer, prng, first_diff_pos));
}

// output a string with 'size' characters, randomness is in [0 26]
std::string generate_random_string(const size_t size, const char randomness)
{
	std::string out;
	out.reserve(size);
	for (size_t i =0; i<size; i++) {
		out.push_back('a'+ std::rand() % (randomness+1));
	};
	return out;
}

TEST(gen_rand_str, a)
{
	ASSERT_EQ( 10UL, generate_random_string(10, 8).size() );
	SPLOG("%s", generate_random_string(10, 8).c_str());
}

void test_encoding(
	const std::string& encoding,
	const std::string& expected_output_encoding,
	size_t number_of_kvs,
	size_t chunk_goal_size)
{
	SPLOG("===================================================");
	SPLOG(" launching test (encoding='%s', # of keys: %lu, chunk size: %lu", encoding.c_str(), number_of_kvs, chunk_goal_size);
	path test_path(make_path("encoding/" + encoding));

	output_stream_params o_stream;
	o_stream.encoding = encoding;
	o_stream.goal_size = chunk_goal_size;

	manifest some_random_kvs;

	std::unique_ptr<kv_sink> sink(o_stream.build(test_path, "encoding_"+encoding, some_random_kvs));

	std::string out_encoding;
	ASSERT_NO_THROW(out_encoding = some_random_kvs.get_encoding());
	ASSERT_EQ(0, expected_output_encoding.compare(out_encoding));

	std::srand(time(0));
	tracked_vector< std::pair<std::string, std::string> > original_data(track_alloc("encoding_test:test_encoding:original_data"));

	char randomness = 8;
	size_t key_size = 10LU;
	size_t value_size = 20LU;

	// bogus
	size_t kv_size = 3LU+key_size+value_size;
	size_t expected_keyvalues_per_chunk = chunk_goal_size / kv_size;
	SPLOG("expected kvs/chunk = %lu", expected_keyvalues_per_chunk);

	for (size_t i=0; i<number_of_kvs; i++) {
		std::string rand_key = generate_random_string(key_size, randomness);
		std::string rand_value = generate_random_string(value_size, randomness);
		original_data.push_back(make_pair(rand_key, rand_value));
		sink->write(rand_key, rand_value);
	}
	sink->close();

	// basic validation of the generated manifest
	size_t actual_num_chunks = 0LU;
	for (auto fi : some_random_kvs) {
		actual_num_chunks++;
		SPLOG("file_info[%s] size: %lu num_records: %lu",
			fi.file.filename().c_str(),
			fi.size,
			fi.num_records
		);
		ASSERT_LE(fi.num_records, expected_keyvalues_per_chunk);
	}
	//ASSERT_EQ( expected_chunks, actual_num_chunks );

	SPLOG("done writing random keys in compressed chunks. Let's read them!");
	manifest_reader mr(some_random_kvs);

	std::string key,value;
	bool read_result = false;
	for (auto kv : original_data) {
		ASSERT_NO_THROW( read_result = mr.read(key, value) );
		ASSERT_TRUE(read_result);
		ASSERT_EQ( 0, key.compare(kv.first) );
		ASSERT_EQ( 0, value.compare(kv.second) );
	}
	ASSERT_FALSE( mr.read(key,value) );

	SPLOG(" SUCCESS ! ");
	SPLOG("===================================================");
	SPLOG(".");
}

TEST_F(encoding, empty)
{
	test_encoding(codec::null, codec::null, 1000*1000, 1000*1000);
}

TEST_F(encoding, force_gzip)
{
	test_encoding("", codec::gzip, 2, 33);
}

TEST_F(encoding, gzip_tiny)
{
	test_encoding(codec::gzip, codec::gzip, 2, 33);
}

TEST_F(encoding, gzip_medium)
{
	test_encoding(codec::gzip, codec::gzip, 1000*1000, 1024*1024);
}

TEST_F(encoding, gzip_big)
{
	test_encoding(codec::gzip, codec::gzip, 4*1000*1000, 64*1024*1024);
}

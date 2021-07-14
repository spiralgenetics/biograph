#include "modules/test/test_utils.h"
#include "gtest/gtest.h"
#include "modules/mapred/path.h"
#include "modules/io/file_io.h"
#include "modules/io/mem_io.h"
#include "modules/io/config.h"

#define CHECK_URL_PARSE(X) ASSERT_EQ(X, path(X).url())
#define CHECK_TYPE(X, X_type) ASSERT_EQ(path(X).type(), X_type)

TEST(path_basic, url_parse)
{
	// check that old paths still work:
	CHECK_URL_PARSE("file:///some/stuff/here");
	CHECK_TYPE("file:///some/stuff/here", path::FILE);
}

TEST(path_basic, append)
{
	path p2 = path("file:///foo").append("bar");
	ASSERT_EQ("file:///foo/bar", p2.url());
}

void test_url(const std::string& url)
{
	path p(url);

	// Write a random string to the file and read it back slowly
	std::string str;
	for (int i = 0; i < 20; i++) {
		str.push_back('A' + rand() % 26);
	}

	p.put(str);
	auto reader = p.read();
	char buf[100];
	size_t r = reader->read(buf, 100);
	ASSERT_EQ(r, (unsigned) 20);
	buf[r] = 0;
	ASSERT_EQ(str, buf);
	ASSERT_EQ(p.get(), str);

	// Write a random string slowly, and read it back as one chunk
	str.clear();
	for (int i = 0; i < 20; i++) {
		str.push_back('A' + rand() % 26);
	}

	auto writer = p.write();
	writer->write(str.data(), 20);
	writer->close();
	std::string rstr = p.get();
	ASSERT_EQ(str, rstr);

	// Test inverted write
	str.clear();
	for (int i = 0; i < 20; i++) {
		str.push_back('A' + rand() % 26);
	}

	mem_io mio("", track_alloc("path_test"));
	mio.write(str.data(), 20);
	auto inverted = p.write_inverted(mio, 20);
	inverted->wait();
	rstr = p.get();

	ASSERT_EQ(str, rstr);
}

TEST(path_basic, rw_file)
{
	try
	{
		auto out_path = std::string("file://") + make_path("path_basic");
		test_url(out_path);
	}
	catch(const io_exception& io)
	{
		printf("GOT AN ERROR: %s\n", io.message().c_str());
		ASSERT_TRUE(0);
	}
}

void test_http_url(const std::string& url, long size, const std::string& filename)
{
	path p(url);
	std::unique_ptr<readable> rd = p.read();
	if(rd == nullptr)
		throw io_exception("Could not get readable from URL");
	file_writer fw(filename);
	io_copy(*rd, fw);
	ASSERT_EQ( size, fsize(filename));
}

template <std::string (*root_path)()>
class path_fixture : public ::testing::Test
{
	public:
		path_fixture() {}

	protected:
		void SetUp() override
		{
			root = path(root_path());
			root.mkdir();
			subdir_1 = root.append("1");
			subdir_1.mkdir();
			subdir_2 = subdir_1.append("2");
			subdir_2.mkdir();
			subdir_3 = subdir_2.append("3");
			subdir_3.mkdir();
			subdir_4 = subdir_3.append("4");
			subdir_4.mkdir();
			file_in_2 = subdir_2.append("foo");
			file_in_2.put("secret stuff");
		}
		void TearDown() override { }

		path root, subdir_1, subdir_2, subdir_3, subdir_4, file_in_2;

		void walk_empty_handler()
		{
			ASSERT_NO_THROW( root.walk( [](const path::walk_params& params) {}));
		}

		void count()
		{
			int directories_count(0), files_count(0);

			root.walk([&](const path::walk_params& params) {
				if (params.state == path::w_dir_enter) {
					directories_count++;
				} else if (params.state == path::w_file) {
					files_count++;
				}
			});

			ASSERT_EQ( files_count,		1 );
			ASSERT_EQ( directories_count,	5 );
		}

		void rm_rf()
		{
			ASSERT_EQ( path::e_directory, root.exists() );
			ASSERT_EQ( path::e_directory, subdir_1.exists() );
			ASSERT_EQ( path::e_directory, subdir_2.exists() );
			ASSERT_EQ( path::e_file, file_in_2.exists() );
			ASSERT_EQ( path::e_directory, subdir_3.exists() );
			ASSERT_EQ( path::e_directory, subdir_4.exists() );

			subdir_1.rmdir(true);

			ASSERT_EQ( path::e_no_exist, subdir_4.exists() );
			ASSERT_EQ( path::e_no_exist, subdir_3.exists() );
			ASSERT_EQ( path::e_no_exist, subdir_2.exists() );
			ASSERT_EQ( path::e_no_exist, subdir_1.exists() );
			ASSERT_EQ( path::e_directory, root.exists() );
		}
};

std::string file_root()	{ return CONF_S(storage_root) + "/file_path_walk"; }

typedef path_fixture<file_root> path_file;

TEST_F(path_file, no_handler)
{
	walk_empty_handler();
}

TEST_F(path_file, count)
{
	count();
}

TEST_F(path_file, rm_rf)
{
	rm_rf();
}

TEST(path, excluded)
{
	ASSERT_FALSE(path(CONF_S(storage_root)).excluded());
	ASSERT_TRUE(path("file:///").excluded());
	ASSERT_TRUE(path("file:///out/spiral/storage").excluded());
	ASSERT_FALSE(path("file:///out/spiral/storage").rmdir(true));
}

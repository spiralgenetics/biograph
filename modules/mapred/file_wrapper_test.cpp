#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <gtest/gtest.h>

#include "modules/io/io.h"
#include "modules/io/log.h"
#include "modules/io/file_wrapper.h"
#include "modules/mapred/path.h"
#include "modules/test/test_utils.h"

TEST(file_wrapper, read_write)
{
	std::string s{"The quick brown fox jumped over the lazy dog.\n"};
	std::string path{make_path("file_wrapper_test")};
	{
		file_wrapper writer(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
		ASSERT_TRUE(writer.is_open());
		ASSERT_NE(writer.get_fd(), -1);
		ssize_t bytes_written = ::write(writer.get_fd(), s.c_str(), s.size());
		ASSERT_EQ(bytes_written, s.size());
		ASSERT_THROW(writer.open(make_path("file_wrapper_test"), O_RDWR), io_exception);
		writer.close();
		ASSERT_FALSE(writer.is_open());
		ASSERT_EQ(writer.get_fd(), -1);
	}
	{
		file_wrapper reader(path, O_RDONLY);
		ASSERT_TRUE(reader.is_open());
		ASSERT_NE(reader.get_fd(), -1);
		std::array<char, 1024> buffer;
		ssize_t bytes_read = ::read(reader.get_fd(), buffer.data(), s.size());
		ASSERT_EQ(bytes_read, s.size());
		reader.close();
		ASSERT_FALSE(reader.is_open());
		ASSERT_EQ(reader.get_fd(), -1);
	}
}


TEST(file_wrapper, error)
{
	file_wrapper some_file;
	ASSERT_FALSE(some_file.is_open());
	ASSERT_EQ(some_file.get_fd(), -1);
	ASSERT_THROW(some_file.open("Does not exist", O_RDONLY), io_exception);
}


TEST(file_wrapper, exists)
{
	path empty_file_path = make_path("empty_file");
	file_wrapper empty_file;
	empty_file.open(empty_file_path.bare_path(), O_RDWR | O_CREAT | O_TRUNC, 0644);
	empty_file.close();
	
	file_wrapper already_exists;
	ASSERT_THROW(already_exists.open(empty_file_path.bare_path(), O_RDWR | O_CREAT | O_EXCL), io_exception);
	ASSERT_NO_THROW(already_exists.open(empty_file_path.bare_path(), O_RDWR | O_CREAT));
}


TEST(file_wrapper, move)
{
	std::string s{"The quick brown fox moved over the lazy dog.\n"};
	std::string path{make_path("file_wrapper_move_test")};

	file_wrapper writer1(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
	ssize_t bytes_written = ::write(writer1.get_fd(), s.c_str(), s.size() / 2);
	ASSERT_EQ(bytes_written, s.size() / 2);
	
	file_wrapper writer2(std::move(writer1));
	ASSERT_FALSE(writer1.is_open());
	ASSERT_EQ(writer1.get_fd(), -1);
	ASSERT_TRUE(writer2.is_open());
	ASSERT_NE(writer2.get_fd(), -1);
	bytes_written = ::write(writer2.get_fd(), s.substr(s.size() / 2).c_str(), s.size() - s.size() / 2);
	ASSERT_EQ(bytes_written, s.size() - s.size() / 2);
	writer2.close();
	
	file_wrapper reader(path, O_RDONLY);
	std::array<char, 1024> buffer;
	ssize_t bytes_read = ::read(reader.get_fd(), buffer.data(), s.size());
	ASSERT_EQ(bytes_read, s.size());
}

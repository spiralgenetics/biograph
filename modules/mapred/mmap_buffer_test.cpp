#include "modules/test/test_utils.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/mmap_vector.h"
#include "modules/io/io.h"
#include "modules/mapred/path.h"
#include <gtest/gtest.h>

TEST(mmap_buffer, basic)
{
	mmap_buffer start(make_path("hello"), 1000);
	strcpy(start.buffer(), "Hello World");
	start.close();
	path p(make_path("hello"));
	std::string r = p.get();
	printf("Size of r = %ld\n", r.size());
	printf("%s\n", r.c_str());
	mmap_buffer restart(make_path("hello"));
	char buf[1000];
	strcpy(buf, restart.buffer());
	printf("%s\n", buf);
	if (strcmp(buf, "Hello World") != 0)
		throw io_exception("Bad buffers");
}

struct some_pod
{
	int m_int;
	double m_double;
	char m_char;
	bool m_bool;
	
	some_pod(int an_int, double a_double, char a_char, bool a_bool)
		: m_int(an_int), m_double(a_double), m_char(a_char), m_bool(a_bool)
	{}
};

bool operator==(const some_pod& lhs, const some_pod& rhs)
{
	return lhs.m_int == rhs.m_int
		&& lhs.m_double == rhs.m_double
		&& lhs.m_char == rhs.m_char
		&& lhs.m_bool == rhs.m_bool;
}

TEST(mmap_vector, basic)
{
	std::string mmap_path = make_path("mmap_vector_mmap");
	
	some_pod element0 = {42, 3.14159, 'p', true};
	some_pod element1 = {101, 2.71828, 'e', true};
	some_pod element2 = {-1, 1.618034, 'f', true};
	some_pod element3 = {65536, 1.414213, '2', true};
	some_pod element4 = {1729, 0.30103, 'L', true};
	
	{
		mmap_vector<some_pod> v(5);
		v.get_buffer().open(mmap_path, v.buffer_size());
		v.push_back(element0);
		v.push_back(element1);
		v.push_back(element2);
		v.push_back(element3);
		v.push_back(element4);
	}
	
	mmap_vector<some_pod> v2(5);
	v2.get_buffer().open(mmap_path);
	ASSERT_EQ(element0, v2[0]);
	ASSERT_EQ(element1, v2[1]);
	ASSERT_EQ(element2, v2[2]);
	ASSERT_EQ(element3, v2[3]);
	ASSERT_EQ(element4, v2[4]);
}

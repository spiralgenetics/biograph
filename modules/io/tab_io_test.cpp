#include "modules/io/tab_io.h"
#include "gtest/gtest.h"

class TypeA
{
public:
	TABD_TRANSFER 
	{
		TFIELD(number);
		TFIELD(string);
	}	

public:
	int number;
	std::string string;	
};

class TypeB
{
public:
	TABD_TRANSFER
	{
		TFIELD(subtype);
		TFIELD(number);
	}

public:
	TypeA subtype;
	unsigned char number;
};

template<class T>
void do_base_test(const T& value) 
{
	T out;
	tabd_deserialize(out, tabd_serialize(value));
	ASSERT_EQ(out, value);
}
	
TEST(tab_io, base_types)
{
	do_base_test<bool>(false);
	do_base_test<bool>(true);
	do_base_test<int>(5);
	do_base_test<int>(5000);
	do_base_test<int>(-5000);
	do_base_test<char>('X');
	do_base_test<char>('\t');
	do_base_test<char>('0');
	do_base_test<unsigned char>('\t');
	do_base_test<unsigned char>(17);
	do_base_test<short>(-4097);
	do_base_test<unsigned short>(1000);
	do_base_test<unsigned int>(0x8fffffff);
	do_base_test<std::string>("Hello");
	do_base_test<std::string>("Hell\to");
	std::string x;
	x.push_back('A');
	x.push_back(0);
	x.push_back('A');
	do_base_test<std::string>(x);
}

TEST(tab_io, exact)
{
	ASSERT_EQ(tabd_serialize(503),"503");
	ASSERT_EQ(tabd_serialize(-503),"-503");
	ASSERT_EQ(tabd_serialize(std::string("Hello")),"Hello");
	ASSERT_EQ(tabd_serialize(std::string("\t")),"\\t");
	ASSERT_EQ(tabd_serialize('X'),"X");
	ASSERT_EQ(tabd_serialize((unsigned char)'X'),"88");
	ASSERT_EQ(tabd_serialize('\t'),"\\t");
}

TEST(tab_io, basic)
{
	TypeA x;
	TypeA y;
	x.number = 10;
	x.string = "Hello";
	std::string tabd = tabd_serialize(x);
	printf("%s\n", tabd.c_str());
	tabd_deserialize(y, tabd);
	ASSERT_EQ(y.number, 10);
	ASSERT_EQ(y.string, "Hello");
}

TEST(tab_io, subtype)
{
	TypeB x;
	TypeB y;
	x.subtype.string = "What Up";
	x.number = 3;
	std::string tabd = tabd_serialize(x);
	printf("%s\n", tabd.c_str());
	tabd_deserialize(y, tabd);
	ASSERT_EQ(y.subtype.string, "What Up");
	ASSERT_EQ(y.number, 3);
}

TEST(tab_io, badtype)
{
	int x;
	char c;
        
	ASSERT_THROW(tabd_deserialize(x, "Q"), io_exception);
	ASSERT_THROW(tabd_deserialize(x, "1\t"), io_exception);
	ASSERT_THROW(tabd_deserialize(x, "\t"), io_exception);
	ASSERT_THROW(tabd_deserialize(c, "\\x"), io_exception);
}



#include "gtest/gtest.h"
#include "modules/io/io.h"
#include "modules/io/log.h"
#include "modules/io/transfer_object.h"
#include "modules/io/json_transfer.h"
#include "modules/io/msgpack_transfer.h"

template<class A>
void round_trip(A& out, const A& in)
{
	std::string json = json_serialize(in);
	printf("As JSON: %s\n", json.c_str());
	A tmp;
	json_deserialize(tmp, json);
	std::string json2 = json_serialize(tmp);
	ASSERT_EQ(json, json2);
	std::string msgpack = msgpack_serialize(tmp);
	printf("As msgpack: ");
	printpack(msgpack);
	printf("\n");
	msgpack_deserialize(out, msgpack);
	std::string msgpack2 = msgpack_serialize(out);
	ASSERT_EQ(msgpack, msgpack2);
}

class TypeA
{
public:
	TRANSFER_OBJECT	
	{
		VERSION(1);
		FIELD(number);
		FIELD(string);
	}	

public:
	int number;
	std::string string;	
};

class TypeB
{
public:
	TRANSFER_OBJECT	
	{
		VERSION(1);
		FIELD(subtype);
		FIELD(number);
	}

public:
	TypeA subtype;
	double number;
};

template<class T>
void do_base_test(const T& in) 
{
	T out;
	round_trip(out, in);
	ASSERT_EQ(out, in);
}
	
TEST(SERIAL, base_types)
{
	do_base_test<bool>(false);
	do_base_test<bool>(true);
	do_base_test<int>(5);
	do_base_test<double>(3.1415);
	do_base_test<std::string>("Hello");
}

TEST(SERIAL, array)
{
	std::vector<int> in;
	std::vector<int> out;

	for(int i = 0; i < 5; i++)
	{
		in.push_back(i*i);
	}

	round_trip(out, in);

	ASSERT_EQ(out.size(), (size_t) 5);
	for(int i = 0; i < 5; i++)
	{
		ASSERT_EQ(i*i, out[i]);
	}
}

TEST(SERIAL, basic)
{
	TypeA in;
	TypeA out;
	in.number = 10;
	in.string = "Hello";
	round_trip(out, in);
	ASSERT_EQ(out.number, 10);
	ASSERT_EQ(out.string, "Hello");
}

TEST(SERIAL, subtype)
{
	TypeB in;
	TypeB out;
	in.subtype.string = "What Up";
	in.number = 3.14;
	round_trip(out, in);
        ASSERT_EQ(out.subtype.string, "What Up");
        ASSERT_EQ(out.number, 3.14);
}

TEST(SERIAL, json_badsyntax)
{
	int x;
	std::string bad = "%junk'";

	ASSERT_THROW(json_deserialize(x, bad), deserialization_error);
}

TEST(SERIAL, json_badtype)
{
        int x;
        std::string bad = json_serialize(3.14);
        
        ASSERT_THROW(json_deserialize(x, bad), deserialization_error);
}

TEST(SERIAL, msgpack_badtype)
{
        int x;
        std::string bad = msgpack_serialize(3.14);
        
        ASSERT_THROW(msgpack_deserialize(x, bad), deserialization_error);
}

TEST(SERIAL, pair)
{
	std::pair<int, double> in;
	std::pair<int, double> out;

	in.first = 2;
	in.second = 3.5;	
	round_trip(out, in);
	ASSERT_EQ(out.first, 2);
	ASSERT_EQ(out.second, 3.5);
}

TEST(SERIAL, map)
{
	std::map<std::string, int> in;
	std::map<std::string, int> out;
	
	in["hello"] = 1;
	in["world"] = 42;
	round_trip(out, in);
	ASSERT_EQ(out.size(), (size_t) 2);
	ASSERT_EQ(out["hello"], 1);
	ASSERT_EQ(out["world"], 42);
}

struct test_flags_orig
{
	TRANSFER_OBJECT 
	{
		VERSION(0);
		FIELD(x);
	}
	int x;	
};

struct test_flags_null
{
	TRANSFER_OBJECT 
	{
		VERSION(0);
		FIELD(x);
		FIELD(y);
	}
	int x;	
	transfer_type_null y;
};

#define TEST_FLAGS_OBJ(name, ...) \
struct test_flags_ ## name \
{ \
	TRANSFER_OBJECT \
	{ \
		VERSION(0); \
		FIELD(x); \
		FIELD(y, ##__VA_ARGS__); \
	} \
	int x; \
	int y; \
}; 

TEST_FLAGS_OBJ(none)
TEST_FLAGS_OBJ(def, 42)
TEST_FLAGS_OBJ(strict, TF_STRICT)
TEST_FLAGS_OBJ(no_def, TF_NO_DEFAULT)
TEST_FLAGS_OBJ(allow_null, TF_ALLOW_NULL)
TEST_FLAGS_OBJ(def_allow_null, 42, TF_ALLOW_NULL)

#define PREP_TEST(test_type) \
	printf("Check test " #test_type "\n"); \
	test_flags_orig tf1; \
	test_flags_ ## test_type tf2; \
	tf1.x = 5; \
	tf2.x = 0; \
	tf2.y = 7; 

#define PREP_NULL(test_type) \
	printf("Check test " #test_type "\n"); \
	test_flags_null tf1; \
	test_flags_ ## test_type tf2; \
	tf1.x = 5; \
	tf2.x = 0; \
	tf2.y = 7; 

#define CHECK_THROW(stype) \
	ASSERT_THROW(stype ## _deserialize(tf2, stype ## _serialize(tf1)), deserialization_error);

#define CHECK_VAL(stype, exp_val) \
	stype ## _deserialize(tf2, stype ## _serialize(tf1)); \
	ASSERT_EQ(tf2.x, 5); \
	ASSERT_EQ(tf2.y, exp_val);

TEST(SERIAL, test_flags)
{
	{ PREP_TEST(none); CHECK_VAL(json, 0); }
	{ PREP_TEST(none); CHECK_VAL(msgpack, 0); }
	{ PREP_TEST(def); CHECK_VAL(json, 42); }
	{ PREP_TEST(def); CHECK_VAL(msgpack, 42); }
	{ PREP_TEST(strict); CHECK_THROW(json); }
	{ PREP_TEST(strict); CHECK_THROW(msgpack); }
	{ PREP_TEST(no_def); CHECK_VAL(json, 7); }
	{ PREP_TEST(no_def); CHECK_VAL(msgpack, 7); }
	{ PREP_TEST(allow_null); CHECK_VAL(json, 0); }
	{ PREP_TEST(allow_null); CHECK_VAL(msgpack, 0); }
	{ PREP_TEST(def_allow_null); CHECK_VAL(json, 42); }
	{ PREP_TEST(def_allow_null); CHECK_VAL(msgpack, 42); }
	{ PREP_NULL(none); CHECK_THROW(json); }
	{ PREP_NULL(none); CHECK_THROW(msgpack); }
	{ PREP_NULL(allow_null); CHECK_VAL(json, 0); }
	{ PREP_NULL(allow_null); CHECK_VAL(msgpack, 0); }
	{ PREP_NULL(def_allow_null); CHECK_VAL(json, 42); }
	{ PREP_NULL(def_allow_null); CHECK_VAL(msgpack, 42); }
}


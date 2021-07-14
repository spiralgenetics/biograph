#include <memory>

#include "modules/test/test_utils.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/task_tree.h"

#include <set> 
#include <boost/filesystem.hpp> 
#include <gtest/gtest.h>

namespace fs = boost::filesystem;

std::set<int> g_result;
std::string g_input;

void check_g_result()
{
	ASSERT_TRUE( g_result.find(4) != g_result.cend() );
	ASSERT_TRUE( g_result.find(3) != g_result.cend() );
	ASSERT_TRUE( g_result.find(2) != g_result.cend() );
	ASSERT_TRUE( g_result.find(1) != g_result.cend() );
	ASSERT_TRUE( g_result.find(0) != g_result.cend() );
}

class task_tree_fixture : public ::testing::Test
{
	public:
		task_tree_fixture() {}

	protected:
		void SetUp() override
		{
			std::string out_dir(make_path("build_ref"));
			test_dir = fs::path(out_dir);
			ref_dir = test_dir / "task_tree";
			task_path = test_dir / "task";
		}
		void TearDown() override
		{
		}

		fs::path test_dir;
		fs::path ref_dir;
		fs::path task_path;

		task_mgr_local tm;

		void check_component_task( std::unique_ptr<task> t )
		{
			int out;
			g_result.clear();
			g_input.clear();
			tm.run_task(out, task_path.native(), std::move(t));
			check_g_result();
		}
};

typedef task_tree_fixture task_tree;

#define make_t(n) \
void test##n( const char* input ) \
{ \
	SPLOG(" in test "#n); \
	SPLOG("input: %s", input); \
	g_result.insert(n); \
	g_input = input; \
} \
typedef leaf_task<test##n> t##n;


make_t(0);
make_t(1);
make_t(2);
make_t(3);
make_t(4);

REGISTER_COMPONENT(t0)
REGISTER_COMPONENT(t1)
REGISTER_COMPONENT(t2)
REGISTER_COMPONENT(t3)
REGISTER_COMPONENT(t4)

TEST_F(task_tree, starts_with)
{
	ASSERT_TRUE( starts_with("foo", "foobar") );
	ASSERT_FALSE( starts_with("bar", "foobar") );
}

TEST_F(task_tree, count_children)
{
	ASSERT_EQ( count_children( {} ), 0 );
	ASSERT_EQ( count_children( { "parallel_11324", "parallel_11324" } ), 1 );
	subtasks_t subtasks { "foo", "bar", "joe", "blow", "serial_12345", "yo", "mama", "serial_12345" };
	ASSERT_EQ( count_children(subtasks), 5 );
}

TEST_F(task_tree, constructor)
{
	std::unique_ptr<serial> s( new serial{ t0(), t1(), t2(), t3(), t4() } );

	for( int i = 0; i<5; i++)
	{
		ASSERT_TRUE( (*s).subtasks[i] == std::string("t"+std::to_string(i)) );
	}

	std::unique_ptr<parallel> p( new parallel{t0(), t1(), t2(), t3(), t4()} );

	for( int i = 0; i<5; i++)
	{
		ASSERT_TRUE( (*p).subtasks[i] == std::string("t"+std::to_string(i)) );
	}

	std::unique_ptr<parallel> p0( new parallel { serial { t1() } } );

	ASSERT_TRUE( (*p0).subtasks.size() == 3UL );
	std::string subtask_0 = (*p0).subtasks[0];
	SPLOG("subtask_0: %s", subtask_0.c_str());
	ASSERT_TRUE( starts_with("serial", subtask_0) );
	ASSERT_TRUE( starts_with("t1", (*p0).subtasks[1]) );
	ASSERT_TRUE( (*p0).subtasks[2] == subtask_0 );

	std::unique_ptr<serial> s0( new serial {
		t0(), parallel {
		t1(),
		t2(),
		serial { t3(), t4() }
		}
	});

	ASSERT_TRUE( starts_with("t0",		(*s0).subtasks[0]) );
	ASSERT_TRUE( starts_with("parallel",	(*s0).subtasks[1]) );
	ASSERT_TRUE( starts_with("t1",		(*s0).subtasks[2]) );
	ASSERT_TRUE( starts_with("t2",		(*s0).subtasks[3]) );
	ASSERT_TRUE( starts_with("serial",	(*s0).subtasks[4]) );
	ASSERT_TRUE( starts_with("t3",		(*s0).subtasks[5]) );
	ASSERT_TRUE( starts_with("t4",		(*s0).subtasks[6]) );
	ASSERT_TRUE( starts_with("serial",	(*s0).subtasks[7]) );
	ASSERT_TRUE( starts_with("parallel",	(*s0).subtasks[8]) );

	ASSERT_TRUE( (*s0).subtasks[4] == (*s0).subtasks[7] );
	ASSERT_TRUE( (*s0).subtasks[1] == (*s0).subtasks[8] );
}

TEST_F(task_tree, serial)
{
	std::unique_ptr<serial> s( new serial { t0(), t1(), t2(), t3(), t4() } );

	std::string input = "task_tree_serial";
	s->input =input;
	check_component_task(std::move(s));
	ASSERT_EQ( input, g_input );
}

TEST_F(task_tree, parallel)
{
	std::unique_ptr<parallel> p(new parallel{ t0(), t1(), t2(), t3(), t4() });

	std::string input = "task_tree_parallel";
	p->input = input;
	check_component_task(std::move(p));
	ASSERT_EQ( input, g_input );
}

TEST_F(task_tree, mixed_basic)
{
	std::unique_ptr<serial> mixed( new serial {parallel()} );

	int out;
	tm.run_task(out, task_path.native(), std::move(mixed));
	ASSERT_EQ( out, 0);
}

TEST_F(task_tree, mixed)
{
	std::unique_ptr<serial> mixed( new serial { parallel { t0(), t1(), t2(), t3(), t4() } } );

	check_component_task(std::move(mixed));
}

TEST_F(task_tree, more_mixed)
{
	std::unique_ptr<parallel> mixed1(
		new parallel {
				serial { t0(), parallel { serial { t1(), t2() } } },
				t3(),
				t4()
		}
	);

	check_component_task(std::move(mixed1));
}


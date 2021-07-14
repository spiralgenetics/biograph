#include "modules/test/test_utils.h"
#include "modules/mapred/task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/taskdb.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"
#include <gtest/gtest.h>

class gen_error_task : public task_impl<gen_error_task>
{
public:
	static std::string s_type() { return "gen_error"; }
	TRANSFER_OBJECT {}
	void run() { throw io_exception("woops!"); }
};

class call_error_task : public task_impl<call_error_task>
{
public:
	call_error_task() : m_first_run(true) {}

	static std::string s_type() { return "call_error"; }
	TRANSFER_OBJECT 
	{
		VERSION(0); 
		FIELD(m_first_run, TF_STRICT); 
		FIELD(m_subtask, TF_STRICT); 
	}

	void run() 
	{
		if (m_first_run) {
			m_subtask = add_subtask(make_unique<gen_error_task>());
			m_first_run = false;
		}
		else {
			set_output(0);
		}
	}
private:
	bool m_first_run;
	subtask_id m_subtask;
};

class factorial_task : public task_impl<factorial_task>
{
public:
	static std::string s_type() { return "factorial"; }

	TRANSFER_OBJECT
	{ 
		VERSION(0);
		FIELD(m_first_run); 
		FIELD(m_input); 
		FIELD(m_subtask); 
	}

	// Default constructor 
	factorial_task() {}

	// Generate initially
	factorial_task(int input) 
		: m_first_run(true)
		, m_input(input) 
	{}

	void run()
	{
		int output = 1;

		printf("Factorial %d: first_run = %d", m_input, m_first_run);
		if (!m_first_run)
			printf(", subtask = '%d'\n",  m_subtask);
		else
			printf("\n");

		if (!m_input)
			set_output(output);
		else {
			if (m_first_run) {
				m_subtask = add_subtask(make_unique<factorial_task>(m_input - 1));
				m_first_run = false;
			}
			else {
				get_output(output, m_subtask);
				output *= m_input;
				set_output(output);
				printf("Output = %d\n", output);
			}
		}
	}

private:
	bool m_first_run;
	int m_input;
	subtask_id m_subtask;
};

REGISTER_TASK(factorial_task);
REGISTER_TASK(gen_error_task);
REGISTER_TASK(call_error_task);

TEST(task, count)
{
	path tmp_path(make_path("task_count"));
	task_mgr_local tm;
	int out;
	tm.run_task(out, tmp_path, make_unique<factorial_task>(5));
	ASSERT_EQ(120, out);
}

TEST(task, cancelworks)
{
	try {
		path tmp_path(make_path("task_cancelworks"));
		task_mgr_local tm;
		tm.run_task(tmp_path, make_unique<call_error_task>());
	}
	catch (const io_exception& e) {
		ASSERT_EQ("woops!", e.message().substr(0, 6));
		return;
	}
	FAIL();
}


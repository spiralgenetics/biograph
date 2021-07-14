#include "modules/test/test_utils.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/taskdb.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"
#include <gtest/gtest.h>

class inner_task : public task_impl<inner_task>
{
public:
	static std::string s_type() { return "inner"; }

	TRANSFER_OBJECT {}

	void run() 
	{
		printf("Running inner task\n");
		int r = 0;
		set_output(r);
	}
};

REGISTER_TASK(inner_task);

class split_task : public task_impl<split_task>
{
public:
	split_task() {}
	split_task(int l) : levels(l) {}

	TRANSFER_OBJECT 
	{
		VERSION(0);
		FIELD(levels, TF_STRICT); 
	}

	static std::string s_type() { return "split"; }

	void run()
	{
		printf("Running split task, level = %d\n", levels);
		if (levels == 0) {
			set_output(levels);
			return;
		}
		levels--;
		if (levels == 0) {
			split_progress(.05, .05);
		}
		else {
			split_progress(.1 * .3, 0.7);
		}
		for (int i = 0; i < 9; i++) {
			add_subtask(make_unique<inner_task>());
		}
	}

private:
	int levels;
};

REGISTER_TASK(split_task);

TEST(progress, showbug)
{
	path tmp_path(make_path("progress"));
	task_mgr_local tm;
	tm.run_task(tmp_path, make_unique<split_task>(2));
}

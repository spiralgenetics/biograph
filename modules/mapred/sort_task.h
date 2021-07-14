#pragma once

#include "modules/mapred/input_stream.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task.h"

#include <vector>

class sort_task : public task_impl<sort_task>
{
public:
	static std::string s_type() { return "sort"; }
	std::string subtype() const override { return reduce; }

	TRANSFER_OBJECT 
	{
		VERSION(0); 
		FIELD(input, TF_STRICT); 
		FIELD(reduce, TF_STRICT); 
		FIELD(reduce_param, TF_STRICT); 
		FIELD(is_summary, TF_STRICT);
		FIELD(goal_size, TF_STRICT); 
		FIELD(update_freq, TF_STRICT); 
		FIELD(max_files, TF_STRICT);
		FIELD(m_sorted, TF_STRICT); 
		FIELD(m_subtasks, TF_STRICT); 
		FIELD(m_round, TF_STRICT);
		FIELD(m_expected_rounds, TF_STRICT);
	}

	void run();
	void prepare(manifest& in, double start);
	void load_results(manifest& out, double end);

	// Presumes this is sorted in parts
	manifest input;
	std::string reduce;
	std::string reduce_param;
	bool is_summary = false;

	size_t goal_size = 64*1024*1024;
	size_t update_freq = 10000;
	size_t max_files = 30;

private:
	manifest m_sorted;
	std::vector<subtask_id> m_subtasks;
	size_t m_round = 0;
	size_t m_expected_rounds;
};

class sorted_reduce_task : public task_impl<sorted_reduce_task>
{
public:
	static std::string s_type() { return "sorted_reduce"; }
	std::string subtype() const override { return reduce; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(reduce, TF_STRICT);
		FIELD(reduce_param, TF_STRICT);
		FIELD(out_sort, TF_STRICT);
		FIELD(presorted, TF_STRICT);
		FIELD(prereduce_goal_size, TF_STRICT);
		FIELD(goal_size, TF_STRICT);
		FIELD(update_freq, TF_STRICT);
		FIELD(m_subtasks, TF_STRICT);
	}

	void run();
	void setup();
	void finish();

	manifest input;  // Presumes this is fully sorted
	std::string reduce;
	std::string reduce_param;
	std::string out_sort;
	bool presorted = false;
	size_t prereduce_goal_size = 8*64*1024*1024;
	size_t goal_size = 64*1024*1024;
	size_t update_freq = 2000;

private:
	std::vector<subtask_id> m_subtasks;
};

#pragma once

#include "modules/mapred/reducer.h"
#include "modules/mapred/input_stream.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task.h"

#include <vector>

class reduce_task : public task_impl<reduce_task>
{
public:
	static std::string s_type() { return "reduce"; }
	std::string subtype() const override { return reduce; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT); 
		FIELD(reduce, TF_STRICT); 
		FIELD(reduce_param, TF_STRICT); 
		FIELD(is_summary, TF_STRICT); 
		FIELD(post_sort, TF_STRICT);
		FIELD(goal_size, TF_STRICT); 
		FIELD(mp_goal_size, TF_STRICT); 
		FIELD(update_freq, TF_STRICT); 
		FIELD(max_files, TF_STRICT);
		FIELD(m_state, TF_STRICT); 
		FIELD(m_subtasks, TF_STRICT); 
	}

	void run();
	void prepare(const manifest& in);
	void prepare_mergepart(const manifest& in);
	void prepare_reduce(const manifest& in);
	void load_results(manifest& out);

	manifest input;

	std::string reduce;
	std::string reduce_param;
	bool is_summary = false;
	std::string post_sort;

	size_t goal_size = 64*1024*1024;
	size_t mp_goal_size = 32*1024*1024;
	size_t update_freq = 10000;
	size_t max_files = 25;

private:
	int m_state = 0;  // 0 = never ran anything, 1 = ran mergepart, 2 = ran reduce
	std::vector<subtask_id> m_subtasks;
};

class reduce_part_task : public task_impl<reduce_part_task>
{
public:
	static std::string s_type() { return "reduce_part"; }
	std::string subtype() const override { return reduce; }

	TRANSFER_OBJECT
	{
		VERSION(0); 
		FIELD(input_stream, TF_STRICT);
		FIELD(output_stream, TF_STRICT);
		FIELD(reduce, TF_STRICT); 
		FIELD(reduce_param, TF_STRICT); 
		FIELD(update_freq, TF_STRICT);
	};

	void run();
	void send_update();

	task_requirements get_requirements() override;

	size_t m_num_proc;
	input_stream_params input_stream;
	output_stream_params output_stream;
	
	std::string reduce;
	std::string reduce_param;

	size_t update_freq;

private:
	void make_reducer();

private:
	std::unique_ptr<reducer> m_reducer;
	std::unique_ptr<kv_source> m_input;
	std::unique_ptr<kv_sink> m_output;
};

#pragma once

#include "modules/mapred/mapper.h"
#include "modules/mapred/input_stream.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task.h"

#include <vector>

class map_task : public task_impl<map_task>
{
public:
	static std::string s_type() { return "map"; }
	std::string subtype() const override { return map; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(map, TF_STRICT);
		FIELD(map_param, TF_STRICT);
		FIELD(stable_sort, TF_STRICT);
		FIELD(input_goal_size, TF_STRICT);
		FIELD(output_goal_size, TF_STRICT);
		FIELD(update_freq, TF_STRICT);
		FIELD(num_partitions, TF_STRICT);
		FIELD(sort, TF_STRICT);
		FIELD(reduce, TF_STRICT);
		FIELD(reduce_param, TF_STRICT);
		FIELD(m_subtasks, TF_STRICT);
		FIELD(is_pipe);
	}

	void run();

	manifest input;

	std::string map;
	std::string map_param;
	bool stable_sort = false;

	size_t input_goal_size = 64*1024*1024;
	size_t output_goal_size = 64*1024*1024;
	size_t update_freq = 1000;
	size_t num_partitions = 1;
	bool is_pipe = false;

	std::string sort;
	std::string reduce;
	std::string reduce_param;

private:
	std::vector<subtask_id> m_subtasks;

	subtask_id make_map_part_task(const input_stream_params& params);
	subtask_id make_map_pipe_task(const input_stream_params & the_input_stream_params);
};

class map_part_task : public task_impl<map_part_task>
{
public:
	static std::string s_type() { return "map_part"; }
	std::string subtype() const override { return map; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input_stream, TF_STRICT);
		FIELD(output_stream, TF_STRICT);
		FIELD(map, TF_STRICT);
		FIELD(map_param, TF_STRICT);
		FIELD(update_freq, TF_STRICT);
	};

	void run();
	void send_update();

	task_requirements get_requirements() override;

	input_stream_params input_stream;
	output_stream_params output_stream;

	std::string map;
	std::string map_param;

	size_t update_freq;

private:
	void make_mapper();

private:
	size_t m_num_proc = 0;
	std::unique_ptr<mapper> m_mapper;
	std::unique_ptr<kv_source> m_input;
	std::unique_ptr<kv_sink> m_output;
};

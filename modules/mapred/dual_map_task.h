#pragma once

#include "modules/mapred/dual_mapper.h"
#include "modules/mapred/input_stream.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task.h"

#include <vector>

class dual_map_task : public task_impl<dual_map_task>
{
public:
	static std::string s_type() { return "dual_map"; }
	std::string subtype() const override { return map; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(map, TF_STRICT);
		FIELD(map_param, TF_STRICT);
		FIELD(input_goal_size, TF_STRICT);
		FIELD(output1_goal_size, TF_STRICT);
		FIELD(output2_goal_size, TF_STRICT);
		FIELD(m_subtasks, TF_STRICT);
	}

	void run();

	manifest input;

	std::string map;
	std::string map_param;

	size_t input_goal_size = 64*1024*1024;
	size_t output1_goal_size = 64*1024*1024;
	size_t output2_goal_size = 64*1024*1024;

private:
	std::vector<subtask_id> m_subtasks;

	subtask_id make_map_part_task(const input_stream_params & the_input_stream_params);
};


class dual_map_part_task : public task_impl<dual_map_part_task>
{
public:
	static std::string s_type() { return "dual_map_part"; }
	std::string subtype() const override { return map; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input_stream, TF_STRICT);
		FIELD(output_stream1, TF_STRICT);
		FIELD(output_stream2, TF_STRICT);
		FIELD(map, TF_STRICT);
		FIELD(map_param, TF_STRICT);
	};

	void run();
	void send_update();

	task_requirements get_requirements() override;

	input_stream_params input_stream;
	output_stream_params output_stream1;
	output_stream_params output_stream2;

	std::string map;
	std::string map_param;

	size_t update_freq;

private:
	void make_mapper();

private:
	size_t m_num_proc = 0;
	std::unique_ptr<dual_mapper> m_mapper;
	std::unique_ptr<kv_source> m_input;
	std::unique_ptr<kv_sink> m_output1;
	std::unique_ptr<kv_sink> m_output2;
};

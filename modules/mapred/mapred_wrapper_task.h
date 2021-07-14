
#pragma once

#include "modules/mapred/pipe_params.h"
#include "modules/mapred/task.h"
#include "modules/mapred/manifest.h"
#include "modules/io/log.h"

class mapred_wrapper_task : public task_impl<mapred_wrapper_task>
{
public:
	static std::string s_type() { return "mapred_wrapper"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(aux_map_inputs, TF_STRICT);
		FIELD(map_params, TF_STRICT);
		FIELD(aux_reduce_inputs, TF_STRICT);
		FIELD(reduce_params, TF_STRICT);
		FIELD(m_state, TF_STRICT);
		FIELD(m_map_task, TF_STRICT);
	};

	void run();

	size_t parts = 4;  // Number of splits for the map

	manifest input;  // The input to be mapped over
	std::vector<manifest> aux_map_inputs;  // Auxilary inputs to map taks
	pipe_params map_params;  // The parameters to the map task

	std::vector<manifest> aux_reduce_inputs;  // Auxilary inputs to map task
	pipe_params reduce_params;  // The parameters to the reduce task

private:
	int m_state = 0;  // 0 = before map, 1 = after map, 2 = after reduce
	std::vector<subtask_id> m_map_task;
	subtask_id m_reduce_task;
};

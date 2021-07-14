#pragma once

#include "modules/mapred/pipe_params.h"
#include "modules/mapred/task.h"
#include "modules/mapred/manifest.h"
#include "modules/io/log.h"

class run_wrapper_task : public task_impl<run_wrapper_task>
{
public:
	static std::string s_type() { return "run_wrapper"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(inputs, TF_STRICT);
		FIELD(num_outputs, TF_STRICT);
		FIELD(params, TF_STRICT);
	};

	void run();

	std::vector<manifest> inputs;
	size_t num_outputs;
	pipe_params params;

	void keep_alive() const;

	task_requirements get_requirements() override
	{
		// Setting cpu_minutes to 60 will allocate one himem worker per task.
		// See mapred/taskdb.cpp
		return task_requirements {
			.profile = "himem",
			.cpu_minutes = 60,
		};
	}
};

#pragma once

#include "modules/mapred/input_stream.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task.h"

#include <vector>

class splitter_task : public task_impl<splitter_task>
{
public:
	static std::string s_type() { return "splitter"; }
	std::string subtype() const override { return splitter; }

	TRANSFER_OBJECT 
	{
		VERSION(0); 
		FIELD(splitter, TF_STRICT); 
		FIELD(map_param, TF_STRICT); 
		FIELD(input, TF_STRICT); 
		FIELD(update_freq, TF_STRICT); 
		FIELD(num_partitions, TF_STRICT);
		FIELD(m_split, TF_STRICT); 
		FIELD(m_subtasks, TF_STRICT); 
	}

	void run();
	void prepare(manifest& in, double start);
	void load_results(manifest& out, double end);

	std::string splitter;
	std::string map_param;

	manifest input;

	size_t update_freq = 10000;
	size_t num_partitions = 1;

private:
	manifest m_split;
	std::vector<subtask_id> m_subtasks;
};

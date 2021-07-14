#pragma once 

#include "modules/mapred/map_task.h"
#include "modules/mapred/reduce_task.h"

class map_reduce_task : public task_impl<map_reduce_task>
{
public:
	static std::string s_type() { return "map_reduce"; }
	std::string subtype() const override { return map + "/" + reduce; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT); 
		FIELD(add_before_reduce, TF_STRICT); 
		FIELD(map, TF_STRICT); 
		FIELD(map_param, TF_STRICT);
		FIELD(sort, TF_STRICT); 
		FIELD(reduce, TF_STRICT); 
		FIELD(reduce_param, TF_STRICT); 
		FIELD(post_sort, TF_STRICT);
		FIELD(is_summary, TF_STRICT); 
		FIELD(use_sort, TF_STRICT);
		FIELD(map_update_freq, TF_STRICT); 
		FIELD(reduce_update_freq, TF_STRICT); 
		FIELD(num_partitions, TF_STRICT); 
		FIELD(prec_map, TF_STRICT);
		FIELD(input_goal_size, TF_STRICT);
		FIELD(temp_goal_size, TF_STRICT); 
		FIELD(mp_goal_size, TF_STRICT); 
		FIELD(output_goal_size, TF_STRICT);
		FIELD(m_state, TF_STRICT); 
		FIELD(m_map_task, TF_STRICT); 
		FIELD(m_reduce_task, TF_STRICT);
	}

	void run();

	manifest input;
	manifest add_before_reduce;

	std::string map;
	std::string map_param;
	std::string sort;
	std::string reduce;
	std::string reduce_param;
	std::string post_sort;

	bool is_summary = false;
	bool use_sort = false;
	size_t map_update_freq = 1000;
	size_t reduce_update_freq = 1000;
	size_t num_partitions = 1;
	double prec_map = 0.5;

	size_t input_goal_size = 64*1024*1024;
	size_t temp_goal_size = 32*1024*1024;
	size_t mp_goal_size = 64*1024*1024;
	size_t output_goal_size = 64*1024*1024;

private:
	int m_state = 0;
	subtask_id m_map_task;
	subtask_id m_reduce_task;
};

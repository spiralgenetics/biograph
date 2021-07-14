
#include "modules/io/transfer_object.h"
#include "modules/mapred/input_stream.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task.h"

// Converts a set of corrected reads small enough to fit into
// memory into a PBWT database, and outputs as a resource
class sort_expand_task : public task_impl<sort_expand_task>
{
public:
	sort_expand_task() 
		: split_size(5000000000l) 
		, m_state(0)
	{}

	static std::string s_type() { return "sort_expand_task"; }

	TRANSFER_OBJECT 
	{ 
		VERSION(0);
		FIELD(split_size, TF_STRICT);
		FIELD(input, TF_STRICT); 
		FIELD(m_state, TF_STRICT);
		FIELD(m_subtasks, TF_STRICT);
	}

	void run();
	void void_progress(double progress) { update_progress(progress); }

	size_t split_size;
	manifest input;

	int m_state;
	std::vector<subtask_id> m_subtasks;
};

// Converts a set of corrected reads small enough to fit into
// memory into a PBWT database, and outputs as a resource
class sort_expand_part_task : public task_impl<sort_expand_part_task>
{
public:
	sort_expand_part_task()
	{}

	static std::string s_type() { return "sort_expand_part_task"; }

        task_requirements get_requirements() override {
		return task_requirements { .profile = "himem", .cpu_minutes = 65 };
        };
	

	TRANSFER_OBJECT 
	{ 
		VERSION(0);
		FIELD(input, TF_STRICT); 
		FIELD(output, TF_STRICT);
	}

	void run();
	void void_progress(double progress) { update_progress(progress); }

	input_stream_params input;
	output_stream_params output;
};


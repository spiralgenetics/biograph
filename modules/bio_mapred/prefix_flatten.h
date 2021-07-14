
#include "modules/io/transfer_object.h"
#include "modules/mapred/input_stream.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task.h"

struct flatten_key {
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(context);
		FIELD(shared);
		FIELD(bits);
	}
	uint8_t context; // How long is the entry
	uint8_t shared;  // How much is shared with previous entry
	uint8_t bits;    // Bit field 1,2,4,8 = A,C,G,T
};

class prefix_flatten_task : public task_impl<prefix_flatten_task>
{
public:
	prefix_flatten_task() 
		: prefix_size(0)
		, m_state(0)
	{}

	static std::string s_type() { return "prefix_flatten_task"; }

	TRANSFER_OBJECT { 
		VERSION(0);
		FIELD(input, TF_STRICT); 
		FIELD(prefix_size, TF_STRICT); 
		FIELD(m_state, TF_STRICT);
		FIELD(m_subtasks, TF_STRICT);
	}

	void run();
	void split();
	void join();
	void void_progress(double progress) { update_progress(progress); }

	manifest input;
	size_t prefix_size;

	int m_state;
	std::vector<subtask_id> m_subtasks;
};

class prefix_flatten_part_task : public task_impl<prefix_flatten_part_task>
{
public:
	prefix_flatten_part_task() 
	{}

	static std::string s_type() { return "prefix_flatten_part_task"; }

	TRANSFER_OBJECT { 
		VERSION(0);
		FIELD(input, TF_STRICT); 
		FIELD(output, TF_STRICT); 
		FIELD(prefix, TF_STRICT); 
	}

	void run();
	void void_progress(double progress) { update_progress(progress); }

	manifest input;
	output_stream_params output;
	dna_sequence prefix;
};


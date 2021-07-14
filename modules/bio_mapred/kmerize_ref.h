
#include "modules/io/transfer_object.h"
#include "modules/mapred/task.h"

struct kmerize_ref_params 
{ 
	kmerize_ref_params() 
	{}

	TRANSFER_OBJECT 
	{ 
		VERSION(0); 
		FIELD(kmer_size);
		FIELD(reference);
	} 
	size_t kmer_size;
	std::string reference;
	void validate() {}  // TODO: Make this better
};

class kmerize_ref_task: public task_impl<kmerize_ref_task>
{
public:
	kmerize_ref_task() : m_state(0)
	{}

	static std::string s_type() { return "kmerize_ref_task"; }

	TRANSFER_OBJECT 
	{ 
		VERSION(0);
		FIELD(params, TF_STRICT);
		FIELD(m_state, TF_STRICT);
		FIELD(m_subtasks, TF_STRICT);
	}

	void run();

	kmerize_ref_params params;

private:
	int m_state;  // 0 - Initial, 1 - Called on supercontigs, 2 - Called sort
	std::vector<subtask_id> m_subtasks;
};

class kmerize_supercontig_task: public task_impl<kmerize_supercontig_task>
{
public:
	kmerize_supercontig_task()
	{}

	static std::string s_type() { return "kmerize_supercontig_task"; }

	TRANSFER_OBJECT 
	{ 
		VERSION(0);
		FIELD(params, TF_STRICT);
		FIELD(the_supercontig, TF_STRICT);
	}

	void run();

	kmerize_ref_params params;
	std::string the_supercontig;
};



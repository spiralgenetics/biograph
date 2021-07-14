
#include "modules/io/transfer_object.h"
#include "modules/mapred/task.h"
#include "modules/mapred/manifest.h"

class compute_coverage_task: public task_impl<compute_coverage_task>
{
public:
	compute_coverage_task()
	{}

	static std::string s_type() { return "compute_coverage_task"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(reference, TF_STRICT);
	}

	void run();
	void void_progress(double progress) { update_progress(progress); }

	manifest input;
	std::string reference;
};

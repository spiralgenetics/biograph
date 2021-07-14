#pragma once

#include "modules/mapred/task.h"
#include "tools/version.h"

class build_reference_task : public task_impl<build_reference_task>
{
public:
	build_reference_task() = default;
	build_reference_task(const std::string& build_output_dir, const std::string& reference_name);

	static std::string s_type() { return "build_reference"; }

	enum STATE
	{
		IMPORT_FASTA,
		MAKE_FLAT,
		MAKE_BWT,
		DONE
	};

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(m_out_dir, TF_STRICT);
		FIELD(m_ref_name, TF_STRICT);
		FIELD(m_state, TF_STRICT);
        FIELD(m_min_n_run, TF_STRICT);
	}

	void run();

	std::string m_out_dir;
	std::string m_ref_name;
	int m_state = IMPORT_FASTA;
  size_t m_min_n_run = 50;
};


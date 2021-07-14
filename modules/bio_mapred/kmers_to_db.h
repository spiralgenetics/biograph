#pragma once

#include "modules/io/transfer_object.h"
#include "modules/mapred/task.h"
#include "modules/mapred/manifest.h"

struct kmers_to_db_params
{
        TRANSFER_OBJECT
        {
                VERSION(0);
                FIELD(ref_name);   // Reference
        }
        std::string ref_name;
        void validate();
};

class kmers_to_db_task: public task_impl<kmers_to_db_task>
{
public:
	kmers_to_db_task()
	{}

	static std::string s_type() { return "kmers_to_db_task"; }

	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(ref_name, TF_STRICT);
	}

	void run();
	void void_progress(double progress) { update_progress(progress); }

	manifest input;
	std::string ref_name;
};


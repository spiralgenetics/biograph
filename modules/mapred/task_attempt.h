#pragma once

#include "modules/mapred/task.h"

struct task_attempt
{
	TRANSFER_OBJECT
	{ 
		VERSION(0);
		FIELD(task_id, TF_STRICT); 
		FIELD(state_counter, TF_STRICT);
		FIELD(attempt, TF_STRICT); 
		FIELD(user, TF_STRICT);
		FIELD(working_path, TF_STRICT); 
		FIELD(type, TF_STRICT); 
		FIELD(state_path, TF_STRICT); 
		FIELD(subtask_outputs, TF_STRICT); 
	}

	std::string task_id; // Global task id
	size_t state_counter;  // Which 'version' of this task is this
	size_t attempt;  // Attempt # for this version
	std::string user; // Which user is this for
	path working_path;  // This is where the task should store any files it creates, unique to an attempt
	std::string type;
	path state_path;
	std::vector<path> subtask_outputs;  // Previous outputs
};

struct subtask_definition
{
	TRANSFER_OBJECT	
	{
		VERSION(0);
		FIELD(id, TF_STRICT); 
		FIELD(type, TF_STRICT);
		FIELD(state_path, TF_STRICT); 
		FIELD(requirements);
		FIELD(subtype);
	}

	subtask_id id;
	std::string type;
	path state_path;
	task_requirements requirements;
	std::string subtype;
};

struct task_attempt_result
{
	enum status
	{
		status_error,
		status_new,
		status_done
	};

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(task_id, TF_STRICT); 
		FIELD(state_counter, TF_STRICT); 
		FIELD(attempt, TF_STRICT); 
		FIELD(result, TF_STRICT); 
		FIELD(cur_part, TF_STRICT); 
		FIELD(future_part, TF_STRICT);
		FIELD(state_path, TF_STRICT); 
		FIELD(output, TF_STRICT); 
		FIELD(error, TF_STRICT); 
		FIELD(duration, TF_STRICT);
		FIELD(subtasks, TF_STRICT);
	}

	std::string task_id; // Global task id
	size_t state_counter; // Which 'version' of this task is this
	size_t attempt;  // Attempt #
	int result; // values defined by enum status
	double cur_part;
	double future_part;
	path state_path;
	path output;
	std::string error;
	size_t duration; // in seconds
	std::vector<subtask_definition> subtasks;
};

task_attempt_result attempt_task(const task_attempt& in);

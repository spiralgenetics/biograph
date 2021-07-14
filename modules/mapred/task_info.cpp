#include "modules/mapred/task_info.h"

task_info::task_info(
	const path& storage, 
	const std::string& user, 
	const std::unique_ptr<task>& task)
	: user(user)
	, subtask_id(0)
	, storage(storage)
	, type(task->type())
	, total_progress(1.0)
	, requirements(task->get_requirements())
	, subtype(task->subtype())
{
	std::string state_postfix = printstring("%06ld_%s_initial_state_%ld", 
		random() % 1000000, _id.c_str(), time(0)
	);
	state_path = storage.append(state_postfix);
	state_path.put(task->get_state());
	init();
}

task_info::task_info(
	const task_info& parent,
	const subtask_definition& subtask,
	int child_id, 
	double tot_progress)
	: _id(parent._id + printstring("_%d", child_id))
	, user(parent.user)
	, pid(parent.pid)
	, parent_id(parent._id)
	, subtask_id(child_id)
	, storage(parent.storage)
	, type(subtask.type)
	, state_path(subtask.state_path)
	, total_progress(tot_progress)
	, requirements(subtask.requirements)
	, subtype(subtask.subtype)
{
	init();
}

std::string task_info::get_subtask_id(int subtask_id)
{
	return _id + printstring("_%d", subtask_id);
}

void task_info::init()
{
	doc_type = "task_info";
	state = ts_ready;
	step = 0;
	error_count = 0;
	subtasks_pending = 0;
	remaining_progress = total_progress;
	progress_children = 0.0;
	prev_progress = 0.0;
	cur_progress = 0.0;
	created = last_update = completed = 0;
}

#include "modules/mapred/task.h"
#include "modules/io/log.h"

std::unique_ptr<task> task::create_task(const std::string& type)
{
	SPLOG_P(LOG_DEBUG, "task::create_task> Creating task %s", type.c_str());
	auto it = task_table().find(type);
	if (it == task_table().end()) {
		return nullptr;
	}
	return it->second();
}

void task::register_type(std::string type, task_builder builder)
{
	if (task_table().count(type) > 0) {
		throw io_exception(printstring("Task type already registered: %s", type.c_str()));
	}
	task_table()[type] = builder;
}

task_requirements task::get_requirements()
{
	return task_requirements {
		.profile = "normal",
		.cpu_minutes = 10,
	};
}

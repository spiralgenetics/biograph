#pragma once

#include "modules/io/transfer_object.h"
#include "modules/mapred/task_info.h"

class path;

class task_mgr
{
public:
	task_mgr(std::shared_ptr<taskdb_iface> db);

	// Adds a new task
	std::string add_job(
		const path& base, 
		std::unique_ptr<task> task,
		const std::string& user
	);

	// -1 = error, 0 = running, 1 = done, throws if task doesn't exist
	int state(const std::string& id) const;
	
	// Cancels a task, task now shows as 'error'
	bool cancel_job(const std::string& id, const std::string& message);
	
	// Removes a job and subtasks, auto cancel if running
	void remove_job(const std::string& id);
	
	// Cancels a task, task now shows as 'error'
	bool resurrect_job(const std::string& id);
	
	// Gets progress for a job/task
	double get_progress(const std::string& id) const;
	
	// Gets the error for a job
	std::string get_error(const std::string& id) const;
	
	// Gets the output for a job
	std::string get_output(const std::string& id) const;

	// Give the full details for a task
	task_info get_task_info(const std::string& id) const;
	
	// Gets next job for user
	bool get_task_for_profile(task_info& ti, const std::string& profile) const;

	summary_result get_summary(const std::string& user) const;

	// Output and deserialization all in one place
	template<class Output>
	void get_output(Output& out, const std::string& id)
	{
		json_deserialize(out, get_output(id));
	}

protected:
	std::shared_ptr<taskdb_iface> m_db;
};

class task_mgr_local
{
public:
	task_mgr_local();

	std::string run_task(const path& base, std::unique_ptr<task> task);

	// This is a nice API for 'local' use
	template<class Output>
	void run_task(Output& out, const path& base, std::unique_ptr<task> task)
	{
		json_deserialize(out, run_task(base, std::move(task)));
	}

private:
	std::shared_ptr<taskdb_iface> m_db;
};

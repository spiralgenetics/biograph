#pragma once

#include "modules/mapred/task.h"
#include "modules/mapred/task_attempt.h"
#include "modules/io/config.h"

class task_runner : public task_context
{
public:
	task_runner(const task_attempt& attempt, task_attempt_result& result)
		: m_attempt(attempt)
		, m_result(result)
	{}

	void run();

	// Used by 'run'
	path generate_name(const std::string& task_unique) const;
	subtask_id add_subtask(std::unique_ptr<task> t) override;
	std::string get_output_string(const subtask_id& id) const override;
	void set_output_string(const std::string& _output) override;
	bool update_progress(double progress) override;  // Default just verifies progress, subclass to send to TM
	void split_progress(double cur_part, double future_part) override;
	path get_root() override { return m_attempt.working_path; }

protected:
	const task_attempt& m_attempt;
	task_attempt_result& m_result;

private:
	double m_progress;
	size_t m_next_subtask;
};

class task_worker;
class update_task_runner : public task_runner
{
public:
	update_task_runner(
		task_worker& tw,
		const task_attempt& attempt,
		task_attempt_result& result
	);

	void run();
	bool update_progress(double progress) override;

private:
	task_worker& m_tw;
	time_t m_last_update = 0;
	int m_update_rate = CONF_T(int, task_update_interval);
	int m_task_timeout = CONF_T(int, task_timeout);
	double m_cur_progress = 0.0;
};

#pragma once

#include <memory>

struct taskdb_iface;
class task_info;
struct task_attempt;
struct task_attempt_result;

class task_worker
{
public:
	task_worker(std::shared_ptr<taskdb_iface> db, int max_errors = 5);

	bool get_attempt_for_profile(task_attempt& ta, const std::string& profile) const;
	bool get_attempt_for_id(task_attempt& ta, const std::string& id) const;

	// Returns false if task is out of date
	bool update_progress(const task_attempt_result& in, double value);

	void apply_results(const task_attempt_result& in);
	void set_max_errors(int max_errors);

private:
	task_attempt convert_attempt(const task_info& ti) const;
	bool handle_error(task_info& ti, const task_attempt_result& in);
	void handle_completion(task_info& ti, const task_attempt_result& in);
	void handle_step(task_info& ti, const task_attempt_result& in);

private:
	std::shared_ptr<taskdb_iface> m_db;
	int m_max_errors;
};

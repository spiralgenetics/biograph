#include "modules/mapred/task_worker.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/task_info.h"
#include "modules/io/log.h"

task_worker::task_worker(std::shared_ptr<taskdb_iface> db, int max_errors)
	: m_db(db)
	, m_max_errors(max_errors)
{}

bool task_worker::get_attempt_for_profile(task_attempt& out, const std::string& profile) const
{
	task_info ti;
	if (!m_db->get_for_profile(ti, profile, false)) {
		return false;
	}
	out = convert_attempt(ti);
	return true;
}

bool task_worker::get_attempt_for_id(task_attempt& out, const std::string& id) const
{
	task_info ti;
	if (!m_db->get_for_worker(ti, id)) {
		return false;
	}
	out = convert_attempt(ti);
	return true;
}

task_attempt task_worker::convert_attempt(const task_info& ti) const
{
	task_attempt ta;
	ta.task_id = ti._id;
	ta.state_counter = ti.step;
	ta.attempt = ti.error_count;
	ta.working_path = ti.storage;
	ta.user = ti.user;
	ta.type = ti.type;
	ta.state_path = ti.state_path;
	ta.subtask_outputs = ti.subtask_outputs;
	return ta;
}

bool task_worker::update_progress(const task_attempt_result& in, double progress)
{
	task_info ti;
	do {
		if (!m_db->get(ti, in.task_id)) {
			return false;  // Task has been removed
		}
		if (ti.step != (int) in.state_counter) {
			return false; // Task has moved on
		}
		if (ti.state != ts_running) {
			return false;  // Not in proper state
		}

		double new_progress = ti.remaining_progress * in.cur_part * progress;
		// SPLOG_P(LOG_DEBUG, "task_worker::update_progress> task type: %s new_progress: %f ti.remaining_progress: %f in.cur_part: %f progress: %f", ti.type.c_str(), new_progress, ti.remaining_progress, in.cur_part, progress);

		if (ti.cur_progress < new_progress) {
			ti.cur_progress = new_progress;
		}
	} while(!m_db->put(ti));
	return true;
}

void task_worker::apply_results(const task_attempt_result& in)
{
	task_info ti;
	// SPLOG_P(LOG_DEBUG, "task_worker::apply_results>");
	do {
		//SPLOG("Getting the current state of result data");
		if (!m_db->get(ti, in.task_id)) {
			// SPLOG("task_worker::apply_results> task has been removed");
			return;
		}
		if (ti.step != (int) in.state_counter) {
			// SPLOG("task_worker::apply_results> task has moved on");
			return;
		}
		if (ti.state != ts_running) {
			// SPLOG("task_worker::apply_results> not in proper state");
			return;
		}
		switch (in.result) {
		case task_attempt_result::status_error:
			if (handle_error(ti, in)) {
				return;
			}
			break;
		case task_attempt_result::status_done:
			handle_completion(ti, in);
			break;
		case task_attempt_result::status_new:
			handle_step(ti, in);
			break;
		}
		//SPLOG("Trying to write new result state");
	} while(!m_db->put(ti));
	// SPLOG("task_worker::apply_results> done");
}

bool task_worker::handle_error(task_info& ti, const task_attempt_result& in)
{
	// SPLOG("task_worker::handle_error>");
	ti.error = in.error;
	ti.error_count++;
	if (ti.error_count > m_max_errors) {
		task_mgr tm(m_db);
		tm.cancel_job(ti.user + "-" + ti.pid, in.error);
		return true;
	}
	ti.state = ts_ready;
	return false;
}

void task_worker::handle_completion(task_info& ti, const task_attempt_result& in)
{
	// SPLOG("task_worker::handle_completion>");
	ti.duration = in.duration;

	// Note: parent is magically updated when state transitions to ts_done
	if (ti.parent_id == "") {
		// It's a master task, run its completion routine
		auto task = task::create_task(ti.type);
		if (!task) {
			return;
		}
		task->load_state(in.state_path.get());
		task->complete(ti, true);
	}

	// Final progress update of 100%
	ti.cur_progress = ti.remaining_progress;

	// Update path
	ti.state_path = in.state_path;
	ti.output_path = in.output;

	// Update state
	ti.state = ts_done;
}

void task_worker::handle_step(task_info& ti, const task_attempt_result& in)
{
	// SPLOG("task_worker::handle_step>");
	// Final progress update of 100%
	ti.cur_progress = ti.remaining_progress * in.cur_part;

	// Compute total progress to be dolled out among children
	ti.progress_children = 0.0;
	if (in.subtasks.size()) {
		ti.progress_children = ti.remaining_progress * (1.0 - in.cur_part - in.future_part);
	}

	// Prepare progress for next step
	ti.prev_progress += ti.cur_progress;
	ti.cur_progress = 0.0;
	ti.remaining_progress = ti.remaining_progress * in.future_part;

	// Add children definitions
	ti.subtasks_definitions.insert(ti.subtasks_definitions.end(), in.subtasks.begin(), in.subtasks.end());
	ti.error_count = 0;
	ti.step++;
	ti.state_path = in.state_path;
	ti.subtasks_pending = in.subtasks.size();
	ti.subtask_outputs.resize(ti.subtask_outputs.size() + in.subtasks.size());
	if (in.subtasks.size()) {
		ti.state = ts_adding_children;
	}
	else {
		ti.state = ts_ready;
	}
}

void task_worker::set_max_errors(int max_errors) {
	m_max_errors = max_errors;
}


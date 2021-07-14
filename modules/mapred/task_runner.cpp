#include "modules/mapred/task_runner.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/task_worker.h"
#include "modules/io/utils.h"
#include "modules/io/log.h"
#include "modules/io/config.h"
#include "modules/io/stopwatch.h"
#include <stdlib.h>

path task_runner::generate_name(const std::string& task_unique) const
{
	int random_prefix = random() % 1000000;
	std::string id = m_attempt.task_id;
	for (size_t i = 0; i < id.size(); i++) {
		if (id[i] == '/') {
			id[i] = '_';
		}
	}
	std::string filename = printstring("%06d_%s_%d_%d_%s_%d",
		random_prefix,
		id.c_str(),
		(int) m_attempt.state_counter,
		(int) m_attempt.attempt,
		task_unique.c_str(),
		(int) time(0)
	);
	return m_attempt.working_path.append(filename);
}

subtask_id task_runner::add_subtask(std::unique_ptr<task> task)
{
	subtask_definition subtask;
	subtask.id = m_next_subtask++;
	subtask.state_path = generate_name(printstring("sdef_%d", subtask.id));
	subtask.type = task->type();
	subtask.subtype = task->subtype();
	subtask.state_path.put(task->get_state());
	subtask.requirements = task->get_requirements();
	m_result.subtasks.push_back(subtask);
	SPLOG_P(LOG_DEBUG, "task_runner::add_subtask> %s", subtask.state_path.url().c_str());
	return subtask.id;
}

std::string task_runner::get_output_string(const subtask_id& id) const
{
	if (id >= m_attempt.subtask_outputs.size()) {
		throw io_exception("Unknown subtask result");
	}
	return m_attempt.subtask_outputs[id].get();
}

void task_runner::set_output_string(const std::string& output)
{
	m_result.output = generate_name("output");
	m_result.output.put(output);
	m_result.result = task_attempt_result::status_done;
}

bool task_runner::update_progress(double progress)
{
	if (progress < 0.0 || progress > 1.0) {
		throw io_exception("Progress out of range");
	}
	if (progress < m_progress) {
          throw io_exception(
              printstring("Progress trying to go backwards from %f to %f", m_progress, progress));
        }
	m_progress = progress;
	return true;
}

void task_runner::split_progress(double cur_part, double future_part)
{
	if (cur_part < 0.0 || future_part < 0.0 || cur_part + future_part > 1.0) {
		throw io_exception("Invalid values passed to split_progress");
	}

	if (m_progress > 0.0 && cur_part < m_result.cur_part) {
		throw io_exception("Trying to reduce cur_part while progress > 0.0");
	}

	m_result.cur_part = cur_part;
	m_result.future_part = future_part;
}

void task_runner::run()
{
	m_result.task_id = m_attempt.task_id;
	m_result.state_counter = m_attempt.state_counter;
	m_result.attempt = m_attempt.attempt;
	m_result.result = task_attempt_result::status_new;
	m_result.cur_part = 1.0;
	m_result.future_part = 0.0;
	m_progress = 0.0;
	m_next_subtask = m_attempt.subtask_outputs.size();

	try {
		auto task = task::create_task(m_attempt.type);
		if (!task) {
			throw io_exception("Unknown task type");
		}
		task->load_state(m_attempt.state_path.get());
		auto milliseconds = stopwatch([&] {
			task->run_task(this);
		});
		m_result.duration = std::chrono::duration_cast<std::chrono::seconds>(milliseconds).count();

		if (m_result.result == task_attempt_result::status_done && m_result.subtasks.size()) {
			throw io_exception("Cannot both return output and make subtasks");
		}
		path out = generate_name("new_state");
		out.put(task->get_state());
		m_result.state_path = out;
	}
	catch (const io_exception& e) {
		SPLOG("Task failed with exception: %s", e.message().c_str());
		m_result.result = task_attempt_result::status_error;
		m_result.error = printstring("%s, state path: %s, task type: %s ",
			e.message().c_str(),
			m_attempt.state_path.filename().c_str(),
			m_attempt.type.c_str()
		);
	}
}

task_attempt_result attempt_task(const task_attempt& in)
{
	task_attempt_result out;
	task_runner runner(in, out);
	runner.run();
	return out;
}

update_task_runner::update_task_runner(
		task_worker& tw,
		const task_attempt& attempt,
		task_attempt_result& result)
	: task_runner(attempt, result)
	, m_tw(tw)
{}

void update_task_runner::run()
{
	try {
		task_runner::run();
	}
	catch (const io_exception& e) {
		SPLOG("Caught an io_exception exception in update_task_runner::run: %s", e.message().c_str());
		m_result.result = task_attempt_result::status_error;
		m_result.error = "An io_exception was thrown during run: ";
		m_result.error += e.message();
	}
	catch (const std::exception& e) {
		SPLOG("Caught a standard exception in update_task_runner::run: %s", e.what());
		m_result.result = task_attempt_result::status_error;
		m_result.error = "A standard exception was thrown during run: ";
		m_result.error += e.what();
	}
	catch (...) {
		SPLOG("Got an unknown exception");
		m_result.result = task_attempt_result::status_error;
		m_result.error = "Unknown c++ exception thrown during run";
	}
}

bool update_task_runner::update_progress(double progress)
{
	time_t now = time(0);
	m_cur_progress = progress;
	if (std::difftime(now, m_last_update) < m_update_rate) {
		return true;
	}

	if (!task_runner::update_progress(progress)) {
		return false;
	}
	
	SPLOG_P(LOG_DEBUG, "update_task_runner::update_progress> Sending update: %f", progress);

	// Do watchdog output
	fprintf(stdout, "U");  // Update
	fflush(stdout);
	fsync(1);

	try {
		auto keep_going = m_tw.update_progress(m_result, progress);
		if (!keep_going) {
			SPLOG("update_task_runner::update_progress> Got terminate signal from taskdb, terminating.");
			std::exit(EXIT_SUCCESS);
		}
		time_t now2 = time(0);
		if (now2 - now > 2) {
			m_update_rate = std::min(m_update_rate * 2, m_task_timeout / 2);
			SPLOG("update_task_runner::update_progress> Backing off update rate. "
				"New update rate is %d sec", m_update_rate);
		}
		else if (m_update_rate > CONF_T(int, task_update_interval)) {
			m_update_rate--;
		}
		else if (m_update_rate < CONF_T(int, task_update_interval)) {
			m_update_rate = CONF_T(int, task_update_interval);
		}
	}
	catch (const io_exception& io) {
		SPLOG("update_task_runner::update_progress> Caught exception %s", io.message().c_str());
		m_update_rate = std::min(m_update_rate * 2, m_task_timeout / 2);
	}
	m_last_update = time(0);

	return true;
}

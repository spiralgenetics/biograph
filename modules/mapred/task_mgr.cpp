#include "base/base.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/task_worker.h"
#include "modules/mapred/taskdb.h"
#include "modules/io/utils.h"

#include <algorithm>
#include <functional>

int summary_to_int(const summary_info& summary)
{
	size_t cpos = summary.job_id.rfind("-");
	return std::stoi(summary.job_id.substr(cpos+1));
}

task_mgr::task_mgr(std::shared_ptr<taskdb_iface> db)
	: m_db(db)
{}

std::string task_mgr::add_job(
	const path& base,
	std::unique_ptr<task> task,
	const std::string& user)
{
	SPLOG_P(LOG_DEBUG, "task_mgr::add_job> start");
	task_info ti(base, user, std::move(task));

	auto summary = get_summary(user);
	std::vector<int> id_list;
	id_list.reserve(summary.size());
	if (summary.empty()) {
		ti.pid = '1';
	} else {
		std::transform(summary.begin(), summary.end(), std::back_inserter(id_list), summary_to_int);
		int maximum = *std::max_element(id_list.begin(), id_list.end(), std::less<int>());
		ti.pid = std::to_string(maximum + 1);
	}
	ti._id = ti.user + "-" + ti.pid;
	m_db->put( ti );
	SPLOG("task_mgr::add_job> generated job ID %s", ti._id.c_str());
	return ti._id;
}

int task_mgr::state(const std::string& id) const
{
	task_info ti;
	if (!m_db->get(ti, id)) {
		return -1;
	}
	if (ti.state == ts_done) {
		return 1;
	}
	if (ti.state == ts_cancelled || ti.state == ts_canceling) {
		return -1;
	}
	return 0;
}

bool task_mgr::cancel_job(const std::string& id, const std::string& message)
{
	SPLOG("task_mgr::cancel_job> cancelling job %s", id.c_str());
	task_info ti;
	do {
		if (!m_db->get(ti, id)) {
			SPLOG("Unable to find job");
			return false;
		}
		if (ti.state == ts_canceling || ti.state == ts_cancelled || ti.state == ts_erasing) {
			SPLOG("Job was in non cancelable stae");
			return true;
		}
		ti.state = ts_canceling;
		ti.error = message;
	} while(!m_db->put(ti));
	return true;
}

void task_mgr::remove_job(const std::string& id)
{
	SPLOG("task_mgr::remove_job> removing job %s", id.c_str());
	task_info ti;
	do {
		if (!m_db->get(ti, id)) {
			return;
		}
		if (ti.state == ts_erasing) {
			return;
		}
		ti.state = ts_erasing;
	} while(!m_db->put(ti));
}

bool task_mgr::resurrect_job(const std::string& id)
{
	SPLOG("task_mgr::resurrect_job> resurrecting job %s", id.c_str());
	task_info ti;
	do {
		if (!m_db->get(ti, id)) {
			return false;
		}
		ti.state = ts_resurrect;
	} while(!m_db->put(ti));
	return true;
}

double task_mgr::get_progress(const std::string& id) const
{
	std::string::size_type cpos = id.rfind("-");
	if (cpos == std::string::npos) {
		throw io_exception("Unable to parse id in get_progress");
	}

	std::string user = id.substr(0, cpos);
	std::string pid = id.substr(cpos+1);
	std::string::size_type upos = pid.find("_");
	if (upos != std::string::npos) {
		pid = pid.substr(0, upos);
	}

	// Get the jobs for user
	summary_key start = { user, pid, id };
	summary_key end = { user, pid, id + "_Z" };
	auto summaries = m_db->find_range("job_info", start, end, 0);
	if (summaries.empty()) {
		return 0.0;
	}

	double progress = summaries[0].progress;
	double goal = summaries[0].progress_goal;
	if (goal == 0.0) {
		return 1.0;
	}
	return progress / goal;
}

task_info task_mgr::get_task_info(const std::string& id) const
{
	task_info ti;
	if (!m_db->get(ti, id)) {
		throw io_exception(printstring("Unknown record in couchdb get for key: %s", id.c_str()));
	}
	return ti;
}

bool task_mgr::get_task_for_profile(task_info& ti, const std::string& profile) const
{
	return m_db->get_for_profile(ti, profile, false);
}

std::string task_mgr::get_error(const std::string& id) const
{
	return get_task_info(id).error;
}

std::string task_mgr::get_output(const std::string& id) const
{
	task_info ti = get_task_info(id);
	if (ti.state != ts_done) {
		throw io_exception(printstring("Trying to get output from unfinished task %s", id.c_str()));
	}
	return ti.output_path.get();
}

summary_result task_mgr::get_summary(const std::string& user) const
{
	summary_key start = { user, "" };
	summary_key end = { user, "Z" };
	return m_db->find_range("job_info", start, end, 0, 2);
}

task_mgr_local::task_mgr_local()
	: m_db(new taskdb())
{}

std::string task_mgr_local::run_task(const path& base, std::unique_ptr<task> task)
{
	task_mgr tm(m_db);
	std::string id = tm.add_job(base, std::move(task), "default");
	while (tm.state(id) == 0) {
		task_info ti;
		bool got_task = tm.get_task_for_profile(ti, "");
		if (!got_task) {
			throw io_exception("Ran out of tasks to do, but super task still isn't done");
		}

		task_worker tw(m_db);

		task_attempt attempt;
		(void)tw.get_attempt_for_id(attempt, ti._id);

		task_attempt_result result = attempt_task(attempt);
		tw.apply_results(result);
	}
	if (tm.state(id) == -1) {
		throw io_exception(tm.get_error(id));
	}

	CHECK(tm.state(id) == 1);
	std::string r = tm.get_output(id);
	tm.remove_job(id);
	return r;
}


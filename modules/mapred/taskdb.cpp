#include "modules/mapred/taskdb.h"
#include "modules/mapred/path.h"
#include "modules/io/log.h"
#include "modules/io/config.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/version.h"
#include "modules/io/pulse.h"

static const size_t TASKDB_VERSION = 4;

std::string taskdb_backup_filename(const std::string& root, const std::string& suffix)
{
	return printstring("%s/taskdb/db.%zu%s",
		root.c_str(), TASKDB_VERSION, suffix.c_str()
	);
}

void operator+=(
	std::map<std::string, task_metric>& lhs,
	const std::map<std::string, task_metric>& rhs)
{
	for (const auto& item : rhs) {
		lhs[item.first].cpu_hours += item.second.cpu_hours;
		lhs[item.first].count += item.second.count;
	}
}

static
void operator+=(
	std::map<std::string, const task_info*>& lhs,
	const std::map<std::string, const task_info*>& rhs)
{
	for (const auto& item : rhs) {
		auto twin = lhs[item.first];
		if (item.second) {
			if ((!twin) ||
				(item.second->state < twin->state) ||
				(item.second->state == twin->state &&
					item.second->last_update < twin->last_update)) {
				lhs[item.first] = item.second;
			}
		}
	}
}

void operator+=(summary_info& lhs, const summary_info& rhs)
{
	for (size_t i = 0; i < task_state_count; i++) {
		lhs.count_states[i] += rhs.count_states[i];
	}
	if (rhs.job_id != "" && rhs.job_id < lhs.job_id) {
		lhs.job_id = rhs.job_id;
	}
	lhs.tasks_by_profile += rhs.tasks_by_profile;
	lhs.metrics_by_profile += rhs.metrics_by_profile;
	lhs.metrics_by_type += rhs.metrics_by_type;
	lhs.progress += rhs.progress;
	if (rhs.progress_goal > lhs.progress_goal) {
		lhs.progress_goal = rhs.progress_goal;
	}
}

void taskdb::map_value(const task_info& ti, summary_key& sk, summary_info& si)
{
	sk = { ti.user, ti.pid, ti._id };
	si.count_states.resize(task_state_count);
	si.count_states[ti.state]++;
	si.job_id = ti._id;
	si.progress = ti.prev_progress + ti.cur_progress;
	si.progress_goal = ti.total_progress;

	task_metric by_profile{double(ti.requirements.cpu_minutes) / 60.0};
	si.metrics_by_profile.insert(std::make_pair(ti.requirements.profile, by_profile));

	if (ti.completed && ti.created) {
		std::string full_type = ti.type;
		if (!ti.subtype.empty()) {
			full_type += "/" + ti.subtype;
		}
		task_metric by_type{ti.duration / 60.0 / 60.0};
		si.metrics_by_type.insert(std::make_pair(full_type, by_type));
	}
	si.tasks_by_profile[ti.requirements.profile] = &ti;
}

void taskdb::persist_global_state()
{
	task_map_t tasks_copy;
	{
		// hold the lock long enough to make an in-memory copy of the current state
		lock_t guard(m_mutex);
		tasks_copy = m_tasks;
	}

	// now write this copied state to more permanent storage
	//
	// first write to a temporary file just in case we run out of disk space or
	// we crash in the middle of writing
	path backup_tmp{taskdb_backup_filename(CONF_S(storage_root), ".tmp")};
	backup_tmp.put(msgpack_serialize(tasks_copy));

	// now move the temporary file to the final backup file
	// this should happen as a single atomic operation
	path backup{taskdb_backup_filename(CONF_S(storage_root))};
	path::move(backup_tmp, backup);
	// SPLOG("taskdb::persist_global_state> Saved global state to %s", backup.url().c_str());
}

void taskdb::restore_global_state()
{
	auto filename = taskdb_backup_filename(CONF_S(storage_root));
	SPLOG("taskdb::restore_global_state> Restoring global state from %s", filename.c_str());
	path backup(filename);
	if (backup.exists() != path::e_file) {
		SPLOG("taskdb::restore_global_state> Taskdb backup not found");
		return;
	}

	// deserialize into a temporary just in case the backup is corrupt or has the wrong version.
	// an exception will be thrown in this case with m_tasks left unmodified.
	task_map_t tasks;
	msgpack_deserialize(tasks, backup.get());

	m_tasks = tasks;
	for (const auto& it : m_tasks) {
		summary_key sk;
		summary_info si;
		map_value(it.second, sk, si);
		if (it.second.state <= ts_pending) {
			m_running_job_info.insert(std::make_pair(sk, si));
		}
		m_job_info.insert(std::make_pair(sk, si));
	}
	SPLOG("taskdb::restore_global_state> Done restoring global state");
}

bool taskdb::do_insert(const std::string& id, task_info& ti)
{
	if (ti._id != "" && ti._id != id) {
		return false;
	}
	if (ti._rev != "") {
		return false;
	}
	ti._id = id;
	ti._rev = "1";
	ti.created = ti.last_update = time(0);
	auto result = m_tasks.insert(std::make_pair(ti._id, ti));
	if (!result.second) {
		return false;
	}
	auto it = result.first;
	summary_key sk;
	summary_info si;
	map_value(it->second, sk, si);
	m_job_info.insert(std::make_pair(sk, si));
	if (ti.state <= ts_pending) {
		m_running_job_info.insert(std::make_pair(sk, si));
	}
	return true;
}

bool taskdb::do_update(task_info& ti)
{
	auto it = m_tasks.find(ti._id);
	if (it == m_tasks.end()) {
		return false;
	}
	if (ti._rev != it->second._rev) {
		return false;
	}
	ti._rev = printstring("%d", atoi(ti._rev.c_str()) + 1);
	it->second = ti;
	it->second.last_update = time(0);
	if (it->second.state > ts_pending) {
		it->second.completed = time(0);
	}
	summary_key sk;
	summary_info si;
	map_value(it->second, sk, si);
	const auto it1 = m_job_info.find(sk);
	const auto it2 = m_running_job_info.find(sk);
	m_job_info.update(it1, si);
	if (it2 == m_running_job_info.end()) {
		if (ti.state <= ts_pending) {
			m_running_job_info.insert(std::make_pair(sk, si));
		}
	}
	else {
		if (ti.state <= ts_pending) {
			m_running_job_info.update(it2, si);
		}
		else {
			m_running_job_info.erase(it2);
		}
	}
	return true;
}

bool taskdb::do_get(task_info& ti, const std::string& id)
{
	auto it = m_tasks.find(id);
	if (it == m_tasks.end()) {
		return false;
	}
	ti = it->second;
	return true;
}

void taskdb::query(
	couch_results<summary_key, summary_info>& result,
	const summary_map_t& map,
	summary_map_t::const_iterator begin,
	summary_map_t::const_iterator end,
	size_t group_level,
	size_t limit) const
{
	couch_row<summary_key, summary_info> row;
	while (limit == 0 || result.rows.size() < limit) {
		if (begin == map.end()) {
			break;
		}
		if (end != map.end() && begin->first >= end->first) {
			break;
		}
		summary_key ks = begin->first;
		while (ks.size() > group_level) {
			ks.pop_back();
		}
		row.key = ks;
		ks.push_back("~");
		auto range_end = map.lower_bound(ks);
		if (range_end == map.end()) {
			range_end = end;
		}
		if (end != map.end() && range_end->first > end->first) {
			range_end = end;
		}
		row.value = map.total(begin, range_end);
		result.rows.push_back(row);
		begin = range_end;
	}
	result.total_rows = m_tasks.size();
	result.offset = 0;
}

bool taskdb::running_work(task_info& ti)
{
	if (m_running_job_info.begin() == m_running_job_info.end()) {
		return false;
	}

	auto si = m_running_job_info.total(
		m_running_job_info.begin(), m_running_job_info.end()
	);

	for (const auto& item : si.tasks_by_profile) {
		auto top = item.second;
		if (top && top->state < (int)ts_ready) {
			ti = *top;
			return true;
		}
	}

	return false;
}

void taskdb::add_children(task_info& ti)
{
	size_t start_subtasks = ti.subtasks_definitions.size() - ti.subtasks_pending;
	double progress_child = ti.progress_children / ti.subtasks_pending;
	for (size_t i = start_subtasks; i < ti.subtasks_definitions.size(); i++) {
		const auto& sd = ti.subtasks_definitions[i];
		task_info subtask(ti, sd, i, progress_child);
		do_insert(subtask._id, subtask);
	}
	ti.state = ts_pending;
	do_update(ti);
}

void taskdb::start_cancel(const std::string& id, const std::string& error)
{
	task_info ti;

	if (!do_get(ti, id)) {
		return;  // Task is already gone, nevermind
	}
	if (ti.state == ts_canceling || ti.state == ts_cancelled || ti.state == ts_erasing) {
		return;  // Task is canceling or erasing already, nevermind
	}
	ti.error = error;
	ti.state = ts_canceling;
	do_update(ti);
}

void taskdb::perform_finalization(task_info& ti)
{
	lock_t guard(m_mutex);
	SPLOG_P(LOG_DEBUG, "taskdb::perform_finalization> Entry");
	try {
		auto task = task::create_task(ti.type);
		if (!task) {
			throw io_exception(printstring("taskdb::perform_finalization> Unknown task type: %s",
				ti.type.c_str()
			));
		}
		else {
			task->load_state(ti.state_path.get());
			task->complete(ti, false);
		}
	}
	catch (const io_exception& io) {
		SPLOG("Failed to finalize: %s", io.message().c_str());
	}
}

void taskdb::propagate_cancel(task_info& ti)
{
	for (size_t i = 0; i < ti.subtasks_definitions.size(); i++) {
		std::string subid = ti.get_subtask_id(i);
		task_info sub_ti;

		if (!do_get(sub_ti, subid)) {
			continue;  // Task is already gone, nevermind
		}
		if (sub_ti.state == ts_canceling || sub_ti.state == ts_cancelled || sub_ti.state == ts_erasing) {
			continue;  // Task is canceling or erasing already, nevermind
		}
		sub_ti.state = ts_canceling;
		do_update(sub_ti);
	}
	ti.state = ts_cancelled;
	do_update(ti);

	if (ti.parent_id == "") {
		perform_finalization(ti);
	}
}

void taskdb::propagate_erase(task_info& ti)
{
	for (size_t i = 0; i < ti.subtasks_definitions.size(); i++) {
		std::string subid = ti.get_subtask_id(i);
		task_info sub_ti;

		if (!do_get(sub_ti, subid)) {
			continue;  // Task is already gone, nevermind
		}
		if (sub_ti.state == ts_erasing) {
			continue;  // Task is erasing already, nevermind
		}
		sub_ti.state = ts_erasing;
		do_update(sub_ti);
	}
	ti.state = ts_cancelled;
	erase(ti);
}

bool taskdb::maybe_work()
{
	task_info ti;
	if (!running_work(ti)) {
		return false;
	}
	switch(ti.state) {
	case ts_erasing:
		propagate_erase(ti);
		break;
	case ts_canceling:
		propagate_cancel(ti);
		break;
	case ts_adding_children:
		add_children(ti);
		break;
	default:
		return false;
	}
	return true;
}

void taskdb::do_pending_work()
{
	while(maybe_work()) {
		// Do nothing, work happens in 'maybe_work'
	}
}

void taskdb::maybe_restart(const task_info& in)
{
	if (in.state != ts_running ||
		time(0) - in.last_update < CONF_T(int, task_timeout)) {
		return;
	}
	task_info ti = in;
	ti.state = ts_ready;
	ti.error = "task timed out";
	ti.error_count++;
	if (ti.error_count > CONF_T(int, task_max_timeouts)) {
		start_cancel(
			ti.user + "-" + ti.pid,
			printstring("Too many errors: Task timed out, state path: %s task type: %s",
				ti.state_path.filename().c_str(),
				ti.type.c_str()
			)
		);
		do_pending_work();
	}
	else {
		do_update(ti);  // If this update fails, no harm will be done
	}
	SPLOG("Time now: %lu", time(0));
	SPLOG("Detected timeout for this task: %s", json_serialize(ti).c_str());
}

void taskdb::update_parent(const task_info& child)
{
	task_info ti;
	if (!do_get(ti, child.parent_id)) {
		return;  // Parent is gone, forget it
	}
	if (ti.state != ts_pending) {
		return;  // Parent isn't waiting for childre
	}
	if (ti.subtask_outputs[child.subtask_id] != path()) {
		return;  // Parent already has a result for this child
	}
	ti.subtask_outputs[child.subtask_id] = child.output_path;
	ti.subtasks_pending--;
	if (ti.subtasks_pending == 0) {
		ti.state = ts_ready;
	}
	do_update(ti);
}

bool taskdb::find_task(task_info& out, const std::string& profile)
{
	// find a task suitable for the specified profile
	// if the profile is an empty string (""), return any available task

	if (m_running_job_info.begin() == m_running_job_info.end()) {
		return false;
	}

	auto job = m_running_job_info.total(
		m_running_job_info.begin(), m_running_job_info.end()
	);

	for (const auto& item : job.tasks_by_profile) {
		if (item.second && (profile.empty() || item.first == profile)) {
			auto top = item.second;
			if (top->state == ts_ready) {
				out = *top;
				return true;
			}
			maybe_restart(*top);
		}
	}

	return false;
}

int taskdb::resurrect_task(task_info& ti)
{
	ti.error_count = 0;
	ti.error = "";
	if (ti.output_path.valid()) {
		ti.state = ts_done;
		do_update(ti);
		return ts_done;
	}
	size_t done_child_count = 0;
	for (size_t i = 0; i < ti.subtasks_definitions.size(); i++) {
		std::string subid = ti.get_subtask_id(i);
		task_info sub_ti;
		do_get(sub_ti, subid);
		int child_state = resurrect_task(sub_ti);
		if (child_state == ts_done) {
			done_child_count++;
		}
	}
	if (ti.subtasks_definitions.size() > 0 && done_child_count < ti.subtasks_definitions.size()) {
		ti.state = ts_pending;
		do_update(ti);
		return ts_pending;
	}
	ti.state = ts_ready;
	do_update(ti);
	return ts_ready;
}

bool taskdb::get_for_profile(task_info& ti, const std::string& profile, bool queued)
{
	lock_t guard(m_mutex);
	if (!find_task(ti, profile)) {
		return false;
	}
	ti.last_update = time(0);
	if (queued) {
		ti.state = ts_queued;
	}
	else {
		ti.state = ts_running;
	}
	do_update(ti);
	return true;
}

bool taskdb::get_for_worker(task_info& ti, const std::string& key)
{
	lock_t guard(m_mutex);
	auto it = m_tasks.find(key);
	if (it == m_tasks.end()) {
		return false;
	}
	ti = it->second;
	ti.last_update = time(0);
	ti.state = ts_running;
	do_update(ti);
	return true;
}

bool taskdb::get(task_info& ti, const std::string& key) const
{
	lock_t guard(m_mutex);
	auto it = m_tasks.find(key);
	if (it == m_tasks.end()) {
		return false;
	}
	ti = it->second;
	ti.current_time = time(0);
	return true;
}

bool taskdb::put(task_info& ti)
{
	lock_t guard(m_mutex);
	if (ti._rev == "") {
		if (!do_insert(ti._id, ti)) {
			return false;
		}
	}
	else {
		if (!do_update(ti)) {
			return false;
		}
		if (ti.state == ts_done && ti.parent_id != "") {
			update_parent(ti);
		}
		if (ti.state == ts_resurrect) {
			resurrect_task(ti);
		}
	}
	do_pending_work();
	return true;
}

bool taskdb::erase(const task_info& doc)
{
	lock_t guard(m_mutex);
	const auto& it = m_tasks.find(doc._id);
	if (it == m_tasks.end()) {
		return false;
	}
	const auto& ti = it->second;
	if (ti._rev != doc._rev) {
		return false;
	}
	summary_key sk;
	summary_info si;
	map_value(ti, sk, si);
	m_job_info.erase(sk);
	m_running_job_info.erase(sk);
	m_tasks.erase(doc._id);
	return true;
}

void taskdb::report(
	const std::string& report_type,
	const summary_key& key,
	const summary_key& startkey,
	const summary_key& endkey,
	size_t group_level,
	size_t limit,
	couch_results<summary_key, summary_info>& result) const
{
	lock_t guard(m_mutex);

	const summary_map_t* map;
	if (report_type == "job_info") {
		map = &m_job_info;
	}
	else {
		map = &m_running_job_info;
	}

	auto begin = map->begin();
	auto end = map->end();
	if (begin == end)
	{
		return;
	}

	if (!key.empty()) {
		begin = map->lower_bound(key);
		end = begin;
		end++;
	}

	if (!startkey.empty()) {
		begin = map->lower_bound(startkey);
	}

	if (!endkey.empty()) {
		end = map->lower_bound(endkey);
	}

	query(result, *map, begin, end, group_level, limit);
}

summary_result taskdb::find_range(
	const std::string& index,
	const summary_key& start,
	const summary_key& end,
	size_t limit,
	int group_level) const
{
	summary_result out;
	summary_key key;
	key.push_back(index);
	couch_results<summary_key, summary_info> result;
	report("job_info", key, start, end, group_level, limit, result);
	for (const auto& row : result.rows) {
		out.push_back(row.value);
	}
	return out;
}

class taskdb_couch : public taskdb_iface
{
public:
	taskdb_couch(const std::string& db_url)
		: m_db(db_url)
	{}

	// taskdb_iface methods
	bool get_for_profile(task_info& ti, const std::string& profile, bool queued) override
	{
		std::string url("by_profile/" + profile);
		if (queued) {
			url += "?queued=true";
		}
		return m_db.get(ti, url);
	}

	bool get_for_worker(task_info& ti, const std::string& id) override
	{
		return m_db.get(ti, "by_worker/" + id);
	}

	bool get(task_info& ti, const std::string& key) const override
	{
		return m_db.get(ti, key);
	}

	bool put(task_info& ti) override
	{
		return m_db.put(ti);
	}

	bool erase(const task_info& ti) override
	{
		return m_db.erase(ti);
	}

	summary_result find_range(
		const std::string& index,
		const summary_key& start,
		const summary_key& end,
		size_t limit,
		int group_level) const override
	{
		return m_db.find_range<summary_info>(index, start, end, limit, group_level);
	}

private:
	couch_server<task_info> m_db;
};

std::shared_ptr<taskdb_iface> new_taskdb_couch()
{
	auto url = make_client_url(
		"taskdb_bind_list",
		"MASTER_PORT_5985_TCP_ADDR",
		"MASTER_PORT_5985_TCP_PORT",
		"/spiral_tasks/"
	);
	auto ptr = new taskdb_couch(url);
	return std::shared_ptr<taskdb_iface>(ptr);
}

class task_profile_handler : public easy_rest_handler
{
public:
	task_profile_handler(taskdb& taskdb, http_request& request)
		: easy_rest_handler(request)
		, m_taskdb(taskdb)
		, m_profile(get_match_result(1))
		, m_queued(request.get_variable("queued", "") == "true")
	{}

	std::string easy_get() override
	{
		task_info ti;
		if (!m_taskdb.get_for_profile(ti, m_profile, m_queued)) {
			return ""; // returns a 204
		}
		return json_serialize(ti);
	}

private:
	taskdb& m_taskdb;
	std::string m_profile;
	bool m_queued;
};

class task_worker_handler : public easy_rest_handler
{
public:
	task_worker_handler(taskdb& taskdb, http_request& request)
		: easy_rest_handler(request)
		, m_taskdb(taskdb)
		, m_id(get_match_result(1))
	{ }

	std::string easy_get() override
	{
		task_info ti;
		if (!m_taskdb.get_for_worker(ti, m_id)) {
			return ""; // returns a 204
		}
		return json_serialize(ti);
	}

private:
	taskdb& m_taskdb;
	std::string m_id;
};

class task_doc_handler : public easy_rest_handler
{
public:
	task_doc_handler(taskdb& taskdb, http_request& request)
		: easy_rest_handler(request)
		, m_taskdb(taskdb)
		, m_id(get_match_result(1))
	{
		try {
			m_rev = request.get_variable("rev");
		}
		catch (const io_exception& io) {
		}
	}

	std::string easy_get() override
	{
		task_info ti;
		if (!m_taskdb.get(ti, m_id)) {
			return ""; // returns a 204
		}
		return json_serialize(ti);
	}

	bool easy_put(const std::string& newvalue) override
	{
		task_info ti;
		json_deserialize(ti, newvalue);
		if (!m_taskdb.put(ti)) {
			throw conflict();
		}
		return true;
	}

	bool easy_del() override
	{
		task_info ti;
		ti._id = m_id;
		ti._rev = m_rev;
		if (!m_taskdb.erase(ti)) {
			throw conflict();
		}
		return true;
	}

private:
	taskdb& m_taskdb;
	std::string m_id;
	std::string m_rev;
};

class task_query_handler : public easy_rest_handler
{
public:
	task_query_handler(taskdb& taskdb, http_request& request)
		: easy_rest_handler(request)
		, m_taskdb(taskdb)
		, type(get_match_result(1))
	{
		key = request.get_variable("key", "");
		startkey = request.get_variable("startkey", "");
		endkey = request.get_variable("endkey", "");
		limit = atoi(request.get_variable("limit", "0").c_str());
		group = (request.get_variable("group", "false") == "true");
		group_level = atoi(request.get_variable("group_level", "0").c_str());
	}

	std::string easy_get() override
	{
		couch_results<summary_key, summary_info> result;
		if (group) {
			group_level = 3;
		}

		summary_key k;
		summary_key ks;
		summary_key ke;
		if (key != "") {
			json_deserialize(k, key);
		}
		if (startkey != "") {
			json_deserialize(ks, startkey);
		}
		if (endkey != "") {
			json_deserialize(ke, endkey);
		}

		m_taskdb.report(type, k, ks, ke, group_level, limit, result);
		return json_serialize(result);
	}

private:
	taskdb& m_taskdb;
	std::string type;
	std::string key;
	std::string startkey;
	std::string endkey;
	size_t limit;
	bool group;
	size_t group_level;
};

void taskdb::register_handlers()
{
	register_handler("/spiral_tasks/by_profile/(.*)", [&](http_request& request){
		return new task_profile_handler(*this, request);
	});
	register_handler("/spiral_tasks/by_worker/(.*)", [&](http_request& request){
		return new task_worker_handler(*this, request);
	});
	register_handler("/spiral_tasks/([^/]*)", [&](http_request& request){
		return new task_doc_handler(*this, request);
	});
	register_handler("/spiral_tasks/view/(.*)", [&](http_request& request){
		return new task_query_handler(*this, request);
	});
}

std::unique_ptr<pulse> g_persister;

void taskdb_start_persister(taskdb& taskdb, std::chrono::seconds period)
{
	g_persister.reset(new pulse(period, [&] (std::cv_status) {
		taskdb.persist_global_state();
	}));
}

void taskdb_stop_persister()
{
	if (g_persister) {
		g_persister->stop();
	}
}

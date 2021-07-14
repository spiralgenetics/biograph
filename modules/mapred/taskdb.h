#pragma once

#include "modules/io/aggregate_map.h"
#include "modules/mapred/task_info.h"
#include "modules/web/couchdb.h"
#include "modules/web/restful.h"

#include <mutex>
#include <memory>
#include <chrono>

typedef std::map<std::string, task_info> task_map_t;
typedef aggregate_map<summary_key, summary_info> summary_map_t;

// for unit testing & garbage_collect
std::string taskdb_backup_filename(const std::string& root, const std::string& suffix = "");

class taskdb : public taskdb_iface
{
public:
	void register_handlers();

	void report(
		const std::string& report_type,
		const summary_key& key, 
		const summary_key& startkey,
		const summary_key& endkey,
		size_t group_level,
		size_t limit,
		couch_results<summary_key, summary_info>& result) const;

	void persist_global_state();
	void restore_global_state();

public:
	// taskdb_iface methods
	bool get_for_profile(task_info& ti, const std::string& profile, bool queued) override;
	bool get_for_worker(task_info& ti, const std::string& id) override;
	bool get(task_info& ti, const std::string& key) const override;
	bool put(task_info& ti) override;
	bool erase(const task_info& ti) override;

	summary_result find_range(
		const std::string& index, 
		const summary_key& start, 
		const summary_key& end, 
		size_t limit = 0, 
		int group_level = 0
	) const override;
	
private:
	void map_value(const task_info& t, summary_key& k, summary_info& i);

	bool do_insert(const std::string& id, task_info& ti);
	bool do_update(task_info& ti);
	bool do_get(task_info& ti, const std::string& id);

	void query(
		couch_results<summary_key, summary_info>& result, 
		const summary_map_t& map,
		summary_map_t::const_iterator begin, 
		summary_map_t::const_iterator end,
		size_t group_level, 
		size_t limit) const;

	bool running_work(task_info& ti);
	void add_children(task_info& ti);
	void start_cancel(const std::string& id, const std::string& error);
	void perform_finalization(task_info& ti);
	void propagate_cancel(task_info& ti);
	void propagate_erase(task_info& ti);
	bool maybe_work();
	void do_pending_work();
	void maybe_restart(const task_info& in);
	void update_parent(const task_info& child);
	bool find_task(task_info& out, const std::string& profile);
	int resurrect_task(task_info& ti);  // Returns state

private:
	// The following variables capture the state of the taskdb.
	// Note, though, that m_job_info and m_running_job_info can be recovered 
	// completely from m_tasks, via map_value.
	task_map_t      m_tasks;
	summary_map_t   m_job_info;
	summary_map_t   m_running_job_info;

	typedef std::recursive_mutex mutex_t;
	typedef std::lock_guard<mutex_t> lock_t;
	mutable mutex_t m_mutex;
};

void taskdb_start_persister(taskdb& taskdb, std::chrono::seconds period);
void taskdb_stop_persister();

#pragma once

#include "base/base.h"
#include "modules/mapred/path.h"
#include "modules/mapred/task_attempt.h"

// Task states
// States 0-2 are 'active' states, in that they require work from the worker
// State 3 is 'semi-active' in that if a worker fails while running the 'running' state must be
// brought back to ready (or cancelled)
// State 8 is also requires handling by task_db
// All other states do not directly require any handling

static const int ts_erasing = 0;         // Task is being erased
static const int ts_canceling = 1;       // Task is being cancelled
static const int ts_adding_children = 2; // Task returned with children, adding subtasks
static const int ts_ready = 3;           // Can be run at anytime
static const int ts_queued = 4;          // Grid Engine has queued this task
static const int ts_running = 5;         // Is being run
static const int ts_pending = 6;         // Is waiting on children
static const int ts_cancelled = 7;       // Cancelled for some reason
static const int ts_done = 8;            // Finished (can be garbage collected)
static const int ts_resurrect = 9;       // Taskdb should resurrect this task, never visibile on reads
static const int ts_invalid = 10;        // Must be the last entry

static const int task_state_count = 11;
static_assert(task_state_count == ts_invalid + 1, "If you change the task_states, be sure and update the count!");

class task_info
{
public:
	task_info() 
		: state(ts_invalid) 
	{}

	// constructor used for creating top-level tasks
	task_info(
		const path& storage, 
		const std::string& user, 
		const std::unique_ptr<task>& task
	);

	// this constructor is used by taskdb when converting subtasks_definitions.
	task_info(
		const task_info& parent, 
		const subtask_definition& subtask,
		int child_id, 
		double total_progress
	);

	std::string get_subtask_id(int subtask_id);

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD_SPECIAL(COUCHDB_RESERVED, _id); 
		FIELD_SPECIAL(COUCHDB_RESERVED, _rev); 
		FIELD(doc_type, TF_STRICT); 
		FIELD(user, TF_STRICT); 
		FIELD(pid, TF_STRICT); 
		FIELD(parent_id, TF_STRICT); 
		FIELD(subtask_id, TF_STRICT); 
		FIELD(storage, TF_STRICT); 
		FIELD(type, TF_STRICT); 
		FIELD(state_path, TF_STRICT); 
		FIELD(state, TF_STRICT); 
		FIELD(step, TF_STRICT);
		FIELD(error_count, TF_STRICT);
		FIELD(output_path, TF_STRICT); 
		FIELD(error, TF_STRICT);
		FIELD(subtasks_pending, TF_STRICT); 
		FIELD(subtask_outputs, TF_STRICT); 
		FIELD(subtasks_definitions, TF_STRICT);
		FIELD(progress_children, TF_STRICT);
		FIELD(total_progress, TF_STRICT); 
		FIELD(remaining_progress, TF_STRICT); 
		FIELD(prev_progress, TF_STRICT); 
		FIELD(cur_progress, TF_STRICT);
		FIELD(created, TF_STRICT); 
		FIELD(last_update, TF_STRICT); 
		FIELD(completed, TF_STRICT);
		FIELD(current_time, TF_STRICT);
		FIELD(requirements, TF_STRICT);
		FIELD(subtype); 
		FIELD(duration); 
	}

	std::string doc_type;
	std::string _id; // Total id, job + subid as string, also id for couchdb
	std::string _rev; // Couchdb revision number
	std::string user; // Needed for keying of cumulative data
	std::string pid; // User specific process id
	std::string parent_id; // Parent id or empty
	int subtask_id;
	path storage; // Where state and data files are written to
	std::string type; // The 'type' of task 
	path state_path; // The path to the task state
	int state; // The state of the task (ts_*)
	int step; // Which 'step' of the task we are on
	int error_count; // How many times we have failed to do this step
	path output_path; // The path of the task output
	std::string error; // The most recent error
	int subtasks_pending; // How many subtasks as still incomplete
	std::vector<path> subtask_outputs; // The output paths of my children (empty for running children)
	std::vector<subtask_definition> subtasks_definitions;  // The definitions of all my children
	double total_progress; // Total weight progress this task will have when complete
	double remaining_progress;  // Total remaing for this and future steps
	double prev_progress; // Weighted progress from prior steps (not including children)
	double cur_progress; // Current progress (weighted to allow immediate addition to prev_progress)
	double progress_children;  // Progress to divide among children (only valid in adding_children state)

	// all the following time stamps are only set by the same unique machine that is hosting a taskdb process
	time_t created; // Creation time
	time_t completed;  // Completion time, initially 0
	time_t last_update; // Time of last update
	// Time on taskdb as of last GET. This is useful to meaningfully compare with (created, completed, last_update),
	// even on a different machine, like the web site or the 'server' processes.
	time_t current_time; 

	task_requirements requirements;
	std::string subtype;
	size_t duration; // in seconds

private:
	void init();
};

struct task_metric
{
	task_metric() = default;

	explicit task_metric(double cpu_hours)
		: cpu_hours(cpu_hours)
		, count(1)
	{}

	TRANSFER_OBJECT
	{ 
		VERSION(0);
		FIELD(cpu_hours); 
		FIELD(count);
	}

	double cpu_hours = 0;
	size_t count = 0;
};

struct summary_info
{
	summary_info() 
		: count_states(task_state_count)
	{}

	TRANSFER_OBJECT
	{ 
		VERSION(0);
		FIELD(count_states); 
		FIELD(job_id); 
		FIELD(progress); 
		FIELD(progress_goal); 
		FIELD(metrics_by_profile);
		FIELD(metrics_by_type);
	}

	std::vector<int> count_states;
	std::string job_id;
	double progress = 0.0;
	double progress_goal = 0.0;
	std::map<std::string, task_metric> metrics_by_profile;
	std::map<std::string, task_metric> metrics_by_type;
	std::map<std::string, const task_info*> tasks_by_profile;
};

typedef std::vector<std::string> summary_key;
typedef std::vector<summary_info> summary_result;

struct taskdb_iface
{
	virtual ~taskdb_iface() {}

	virtual bool get_for_profile(task_info& ti, const std::string& profile, bool queued) = 0;
	virtual bool get_for_worker(task_info& ti, const std::string& id) = 0;
	virtual bool get(task_info& ti, const std::string& key) const = 0;
	virtual bool put(task_info& ti) = 0;
	virtual bool erase(const task_info& ti) = 0;

	virtual summary_result find_range(
		const std::string& index, 
		const summary_key& start, 
		const summary_key& end, 
		size_t limit = 0, 
		int group_level = 0
	) const = 0;
};

std::shared_ptr<taskdb_iface> new_taskdb_couch();

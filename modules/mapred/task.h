#pragma once

#include "modules/mapred/path.h"
#include "modules/io/transfer_object.h"
#include "modules/io/make_unique.h"

#include <memory>

// Passed to the run method of a task
class task;

typedef unsigned int subtask_id;

struct task_requirements
{
	TRANSFER_OBJECT	
	{
		VERSION(0);
		FIELD(profile, TF_STRICT); 
		FIELD(cpu_minutes, TF_STRICT);
	}

	std::string profile;
	size_t cpu_minutes;
};

class task_context
{
public:
	virtual ~task_context() {}
	virtual subtask_id add_subtask(std::unique_ptr<task> t) = 0;
	virtual std::string get_output_string(const subtask_id& id) const = 0;
	virtual void set_output_string(const std::string& output) = 0;
	// Splits the progress between this task (cur part), next child procs (1 - cp - fp), 
	// and everything else (future_part)
	virtual void split_progress(double cur_part, double future_part) = 0;
	// Update the 'current' progress
	virtual bool update_progress(double progress) = 0;
	virtual path get_root() = 0;
};

class task_manager;
class task_info;

class task
{
	friend class task_manager;
public:
	virtual ~task() {}
	static std::unique_ptr<task> create_task(const std::string& type);

	virtual std::string type() const = 0;
	virtual std::string subtype() const { return ""; }
	virtual std::string get_state() const = 0;
	virtual void load_state(const std::string& state) = 0;
	virtual void run_task(task_context* context) = 0;
	virtual void complete(const task_info& ti, bool success) {}  // Runs when task is finalized for any reason (success or failure)

	typedef std::unique_ptr<task> (*task_builder)();
	static void register_type(std::string type, task_builder builder);

	virtual task_requirements get_requirements();

private:
	static std::map<std::string, task_builder>& task_table() 
	{ 
		static std::map<std::string, task_builder> map;
		return map;
	}
};

template<class Derived>
class task_impl : public task
{
public:
	// TODO: Magic to auto register etc
	std::string type() const override { return Derived::s_type(); }

	std::string get_state() const override
	{
		return json_serialize(*static_cast<const Derived*>(this)); 
	}

	void load_state(const std::string& state) override
	{
		json_deserialize(*static_cast<Derived*>(this), state); 
	}

	void run_task(task_context* context) override
	{
		m_task_context = context;
		static_cast<Derived*>(this)->run();
	}

protected:
	template<class Output>
	void set_output(const Output& out)
	{
		m_task_context->set_output_string(json_serialize(out));
	}

	std::string get_output_string(const subtask_id& task_id) 
	{
		return m_task_context->get_output_string(task_id);
	}

	template<class Output>
	void get_output(Output& out, const subtask_id& task_id)
	{
		json_deserialize(out, m_task_context->get_output_string(task_id));
	}

	subtask_id add_subtask(std::unique_ptr<task> t) 
	{
		return m_task_context->add_subtask(std::move(t)); 
	}

	void split_progress(double cur_part, double future_part) 
	{ 
		return m_task_context->split_progress(cur_part, future_part); 
	}

	// progress is between 0.0 and 1.0
	bool update_progress(double progress) const 
	{
		return m_task_context->update_progress(progress); 
	}

	path get_root() { return m_task_context->get_root(); }

private:
	task_context* m_task_context = nullptr;
};

template<class Derived>
class task_register_helper
{
public:
	static std::unique_ptr<task> build_task() { return make_unique<Derived>(); }
	task_register_helper()
	{
		task::register_type(Derived::s_type(), build_task); 
	}
};

#define REGISTER_TASK(cls) \
	static task_register_helper<cls> _register_task ## cls;

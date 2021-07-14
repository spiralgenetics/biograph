
#include "modules/mapred/reduce_task.h"
#include "modules/io/utils.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/sorter.h"

REGISTER_TASK(reduce_task);
REGISTER_TASK(reduce_part_task);

const size_t k_part_multiply = 14;
const double k_goal_fullness = .7;

void reduce_task::run()
{
	if (m_state == 0) {
		if (input.get_num_records() == 0) { m_state = 2; return; }
		prepare(input);
	}
	else if (m_state == 1) {
		split_progress(0.02, 0.5);  // Default until overridden by prepare
		manifest temp(input.get_sort(), 1);  // Paritions will adjust during load
		load_results(temp);
		prepare(temp);
	}
	else {
		manifest out(post_sort);
		load_results(out);
		set_output(out);
	}
}

void reduce_task::prepare(const manifest& in)
{
	size_t worst_files = in.max_files();
	if (worst_files <= max_files) {
		split_progress(0.02, 0.02);
		prepare_reduce(in);
		m_state = 2;
	}
	else {
		split_progress(0.02, 0.85);
		prepare_mergepart(in);
		m_state = 1;
	}
}

void reduce_task::prepare_mergepart(const manifest& in)
{
	std::vector<input_stream_params> inputs;
	size_t total_goal_size = size_t(double(mp_goal_size * k_part_multiply) * k_goal_fullness);
	in.split_mergepart(inputs, total_goal_size, max_files);

	for (size_t i = 0; i < inputs.size(); i++) {
		double prec_prog = double(i)  / double(inputs.size());
		update_progress(m_state == 0 ? prec_prog : .5 + .5 * prec_prog);
		auto task = make_unique<reduce_part_task>();
		task->input_stream = inputs[i];
		task->output_stream.goal_size = mp_goal_size;
		task->output_stream.num_partitions = in.get_num_partitions() * k_part_multiply;
		task->output_stream.presorted = true;
		task->output_stream.sort = in.get_sort();
		task->update_freq = update_freq;
		if (is_summary) {
			task->reduce = reduce;
			task->reduce_param = reduce_param;
		}
		else {
			task->reduce = "identity";
		}
		m_subtasks.push_back(add_subtask(std::move(task)));
	}
}

void reduce_task::prepare_reduce(const manifest& in)
{
	std::vector<input_stream_params> inputs;
	in.split_by_partition(inputs);

	for (size_t i = 0; i < inputs.size(); i++) {
		double prec_prog = double(i)  / double(inputs.size());
		update_progress(m_state == 0 ? prec_prog : .5 + .5 * prec_prog);
		auto task = make_unique<reduce_part_task>();
		task->input_stream = inputs[i];
		task->output_stream.goal_size = goal_size;
		task->output_stream.sort = post_sort;
		task->reduce = reduce;
		task->reduce_param = reduce_param;
		task->update_freq = update_freq;
		m_subtasks.push_back(add_subtask(std::move(task)));
	}
}

void reduce_task::load_results(manifest& out)
{
  std::unique_ptr<reducer> reduceit = reducer_registry::get(reduce, reduce_param);
	if (reduceit == nullptr) {
		throw io_exception("Unknown reducer: " + reduce);
	}

	std::string meta = reduceit->get_meta();
	for (size_t i = 0; i < m_subtasks.size(); i++) {
		manifest subout;
		get_output(subout, m_subtasks[i]);
		meta = reduceit->combine_meta(meta, subout.get_meta());
		out.add(subout);
		double prec_prog = double(i)  / double(m_subtasks.size());
		update_progress(m_state == 1 ? .5 * prec_prog : prec_prog);
	}
	out.set_meta(meta);
	m_subtasks.clear();
}

void reduce_part_task::make_reducer()
{
	if (!m_reducer) {
		m_reducer = reducer_registry::get(reduce, reduce_param);
		if (!m_reducer) {
			throw io_exception(printstring("Unknown reducer: %s", reduce.c_str()));
		}
	}
}

task_requirements reduce_part_task::get_requirements()
{
	make_reducer();
	return m_reducer->get_requirements();
}

void reduce_part_task::send_update()
{
	double val = double(m_num_proc) / double(input_stream.num_records);
	update_progress(val);
}

void reduce_part_task::run()
{
	SPLOG_P(LOG_DEBUG, "reduce_part_task::run> Making sorter and reducer");
	std::unique_ptr<sorter> sorter_ = sorter_registry::get(input_stream.sort, "");
	if (!sorter_) {
		throw io_exception("Unknown sorter: " + input_stream.sort);
	}

	make_reducer();

	manifest mout;
	// SPLOG("reduce_part_task::run> Building Input Stream");
	m_input = input_stream.build();
	// SPLOG("reduce_part_task::run> Building Output Stream");
	m_output = output_stream.build(get_root(), "reduce", mout);

	// Do the actual work
	std::string prev_key;
	std::string key;
	std::string value;
	bool first = true;
	m_num_proc = 0;
	// SPLOG("reduce_part_task::run> Entering main loop");

	m_reducer->set_watchdog(std::bind(&reduce_part_task::send_update, this));

	while (m_input->read(key, value)) {
		//SPLOG("Key = %s", key.c_str());
		if (first) {
			m_reducer->start(key, *m_output);
			prev_key = key;
			first = false;
		}
		else  {
			int res = sorter_->compare(prev_key, key);
			if (res > 0) {
				throw io_exception("Key inversion '" + prev_key + "' >= '" + key + "'");
			}
			if (res == -2) {
				m_reducer->end(*m_output);
				m_reducer->start(key, *m_output);
			}
			prev_key = key;
		}
		m_reducer->add_value(key, value, *m_output);
		m_num_proc++;
		send_update();
	}
	if (!first) {
		m_reducer->end(*m_output);
	}
	m_reducer->finalize(*m_output);

	SPLOG_P(LOG_DEBUG, "reduce_part_task::run> Closing output, %ld records processed.", m_num_proc);
	m_output->close();

	// SPLOG("reduce_part_task::run> Writing manifest");
	mout.set_meta(m_reducer->get_meta());
	set_output(mout);
	SPLOG_P(LOG_DEBUG, "reduce_part_task::run> %ld reduce records found.", mout.get_num_records());

	// SPLOG("reduce_part_task::run> Done");
}

#include "modules/mapred/reduce_task.h"
#include "modules/mapred/sort_task.h"
#include "modules/mapred/sorter.h"
#include "modules/io/utils.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"

REGISTER_TASK(sort_task);
REGISTER_TASK(sorted_reduce_task);

void sort_task::run()
{
	split_progress(0.01, 0.5); // Temp guess to make local progress reasonable
	if (m_subtasks.size() == 0) {
		SPLOG("sort_task::run> stage 1: creating reduce_part_tasks ");
		m_round = 1;
		m_expected_rounds = 0;
		size_t cnt = input.get_size();
		while (cnt > 0) {
			cnt /= max_files;
			m_expected_rounds++;
		}
		prepare(input, 0.0);
		return;
	}
	SPLOG("sort_task::run> stage 2: gathering the results of reduce_part_tasks");
	manifest old_out(input.get_sort(), 1);
	load_results(old_out, 0.5);
	prepare(old_out, 0.5);
}

void sort_task::prepare(manifest& in, double start)
{
	std::vector<input_stream_params> inputs;
	m_sorted = manifest(input.get_sort(), 1);
	in.split_sort(m_sorted, inputs, max_files, is_summary);
	m_sorted.merge_tags(in);
	//SPLOG("Sort: inputs.size() = %d, m_sorted.num_records = %d", (int) inputs.size(), (int) m_sorted.get_num_records());
	if (inputs.size() == 0) {
		set_output(m_sorted);
		return;
	}
	double rest = 0.3;
	if (m_round < m_expected_rounds) {
		rest = double(m_expected_rounds - m_round) / double(m_expected_rounds - m_round + 1);
	}
	split_progress(0.01, rest);

	for (size_t i = 0; i < inputs.size(); i++) {
		std::unique_ptr<reduce_part_task> t = make_unique<reduce_part_task>();
		t->input_stream = inputs[i];
		t->output_stream.goal_size = goal_size;
		t->output_stream.num_partitions = 1;
		t->output_stream.presorted = true;
		t->output_stream.sort = in.get_sort();
		t->update_freq = update_freq;
		if (is_summary) {
			t->reduce = reduce;
			t->reduce_param = reduce_param;
		}
		else {
			t->reduce = "identity";
			t->output_stream.allow_split = true;
		}
		m_subtasks.push_back(add_subtask(std::move(t)));
		update_progress(start + (1.0 - start)*i/input.get_size());
	}
}

void sort_task::load_results(manifest& out, double end)
{
	out.add(m_sorted);
	for (size_t i = 0; i < m_subtasks.size(); i++) {
		manifest subout;
		get_output(subout, m_subtasks[i]);
		out.add(subout);
		update_progress(end*i/m_subtasks.size());
	}
	m_subtasks.clear();
	m_round++;
}

void sorted_reduce_task::run()
{
	if (input.get_size() == 0) {
		set_output(input);
		return;
	}
	if (m_subtasks.size() == 0) {
		setup();
	}
	else {
		finish();
	}
}

void sorted_reduce_task::setup()
{
	std::unique_ptr<sorter> sorter_ = sorter_registry::get(input.get_sort(), "");
	std::vector<input_stream_params> inputs;
	input.split_sort_reduce(inputs, prereduce_goal_size);
	split_progress(0.01, 0.01);

	for (size_t i = 0; i < inputs.size(); i++) {
		std::unique_ptr<reduce_part_task> t= make_unique<reduce_part_task>();
		t->input_stream = inputs[i];
		t->output_stream.goal_size = goal_size;
		t->output_stream.num_partitions = 1;
		t->output_stream.presorted = presorted;
		if (presorted) {
			t->output_stream.begin_on = t->input_stream.begin_on;
		}
		t->output_stream.sort = out_sort;
		// Slightly hacky
		t->input_stream.begin_on = sorter_->bump_back(t->input_stream.begin_on);

		t->update_freq = update_freq;
		t->reduce = reduce;
		t->reduce_param = reduce_param;
		m_subtasks.push_back(add_subtask(std::move(t)));
	}
}

void sorted_reduce_task::finish()
{
	manifest out(out_sort);
	for (size_t i = 0; i < m_subtasks.size(); i++) {
		manifest subout;
		get_output(subout, m_subtasks[i]);
		out.add(subout);
		update_progress(double(i) / double(m_subtasks.size()));
	}
	set_output(out);
}

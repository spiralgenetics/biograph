#include "modules/mapred/map_task.h"
#include "modules/mapred/sorter.h"
#include "modules/io/utils.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"

REGISTER_TASK(map_task);
REGISTER_TASK(map_part_task);

void map_task::run()
{
	if (m_subtasks.size() == 0 && input.get_num_records() != 0) {
		SPLOG_P(LOG_DEBUG, "map_task::run> Stage 1");
		split_progress(0.02, 0.05);
		std::vector<input_stream_params> inputs;
		input.split_by_goal_size(inputs, input_goal_size);

		for (size_t i = 0; i < inputs.size(); i++) {
            if (is_pipe) {
                    m_subtasks.push_back(make_map_pipe_task(inputs[i]));
            }
            else {
					m_subtasks.push_back(make_map_part_task(inputs[i]));
            }
			update_progress(double(i) / double(inputs.size()));
		}
		return;
	}

	SPLOG_P(LOG_DEBUG, "map_task::run> Stage 2");
	std::string out_sort = sort;
	if (stable_sort && input.get_sort() != "" && sort == "") {
		out_sort = input.get_sort();
	}
	manifest out(out_sort, num_partitions);

	for (size_t i = 0; i < m_subtasks.size(); i++) {
		manifest subout;
		get_output(subout, m_subtasks[i]);
		out.add(subout);
		update_progress(.85 * double(i) / double(m_subtasks.size()));
	}
	set_output(out);
}


subtask_id map_task::make_map_part_task(const input_stream_params & the_input_stream_params)
{
	auto task = make_unique<map_part_task>();

	task->input_stream = the_input_stream_params;
	task->output_stream.goal_size = output_goal_size;
	task->output_stream.num_partitions = num_partitions;
	if (stable_sort && input.get_sort() != "" && sort == "") {
		task->output_stream.presorted = true;
		task->output_stream.sort = input.get_sort();
	}
	else {
		task->output_stream.sort = sort;
	}

	task->output_stream.reduce = reduce;
	task->output_stream.reduce_param = reduce_param;
	task->map = map;
	task->map_param = map_param;
	task->update_freq = update_freq;

	return add_subtask(std::move(task));
}

task_requirements map_part_task::get_requirements()
{
	make_mapper();
	return m_mapper->get_requirements();
}

void map_part_task::make_mapper()
{
	if (!m_mapper) {
		m_mapper = mapper_registry::get(map, map_param);
		if (!m_mapper) {
			throw io_exception(printstring("Unknown mapper: %s", map.c_str()));
		}
	}
}

void map_part_task::send_update()
{
	double val = double(m_num_proc) / double(input_stream.num_records);
	if (val > 1.0) {
		val = 1.0;
	}
	update_progress(val);
}

void map_part_task::run()
{
	manifest mout;
	SPLOG_P(LOG_DEBUG, "map_part_task::run> Making mapper %s", map.c_str());
	m_num_proc = 0;
	make_mapper();
	m_mapper->set_watchdog(std::bind(&map_part_task::send_update, this));
	m_mapper->setup();

	SPLOG_P(LOG_DEBUG, "map_part_task::run> Building Input Stream");
	m_input = input_stream.build();
	SPLOG_P(LOG_DEBUG, "map_part_task::run> Building Output Stream");
	m_output = output_stream.build(get_root(), "map",  mout);

	// Do the actual work
	std::string key;
	std::string value;
	SPLOG_P(LOG_DEBUG, "map_part_task::run> Entering main loop");
	while (m_input->read(key, value)) {
		m_mapper->map(key, value, *m_output);
		m_num_proc++;
		send_update();
	}

	m_mapper->install_metadata(mout.metadata());

	SPLOG_P(LOG_DEBUG, "map_part_task::run> Closing output");
	m_output->close();

	SPLOG_P(LOG_DEBUG, "map_part_task::run> Writing manifest");
	set_output(mout);

	SPLOG_P(LOG_DEBUG, "map_part_task::run> Done");
}


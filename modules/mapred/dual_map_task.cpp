#include "modules/mapred/dual_map_task.h"
#include "modules/io/utils.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"

REGISTER_TASK(dual_map_task);
REGISTER_TASK(dual_map_part_task);

void dual_map_task::run()
{
	if (m_subtasks.size() == 0 && input.get_num_records() != 0) { // Stage 1
		split_progress(0.02, 0.05);
		std::vector<input_stream_params> inputs;
		input.split_by_goal_size(inputs, input_goal_size);

		for (size_t i = 0; i < inputs.size(); i++) {
			m_subtasks.push_back(make_map_part_task(inputs[i]));
			update_progress(double(i) / double(inputs.size()));
		}
		return;
	}

	// State 2
	std::vector<manifest> out(2);
	for (size_t i = 0; i < m_subtasks.size(); i++) {
		std::vector<manifest> subout;
		get_output(subout, m_subtasks[i]);
		out[0].add(subout[0]);
		out[1].add(subout[1]);
		update_progress(.85 * double(i) / double(m_subtasks.size()));
	}
	set_output(out);
}


subtask_id dual_map_task::make_map_part_task(const input_stream_params & the_input_stream_params)
{
	auto task = make_unique<dual_map_part_task>();

	task->input_stream = the_input_stream_params;
	task->output_stream1.goal_size = output1_goal_size;
	task->output_stream1.num_partitions = 1;
	task->output_stream2.goal_size = output2_goal_size;
	task->output_stream2.num_partitions = 1;
	task->map = map;
	task->map_param = map_param;

	return add_subtask(std::move(task));
}

void dual_map_part_task::make_mapper()
{
	if (!m_mapper) {
		m_mapper = dual_mapper_registry::get(map, map_param);
		if (!m_mapper) {
			throw io_exception(printstring("Unknown mapper: %s", map.c_str()));
		}
	}
}

task_requirements dual_map_part_task::get_requirements()
{
	make_mapper();
	return m_mapper->get_requirements();
}

void dual_map_part_task::send_update()
{
	double val = double(m_num_proc) / double(input_stream.num_records);
	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::send_update> Sending update: %f", val);

	if (val > 1.0) {
		val = 1.0;
	}
	update_progress(val);
}

void dual_map_part_task::run()
{
	std::vector<manifest> out(2);
	SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Making mapper %s", map.c_str());
	m_num_proc = 0;
	make_mapper();
	m_mapper->set_watchdog(std::bind(&dual_map_part_task::send_update, this));
	m_mapper->setup();

	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Building Input Stream");
	m_input = input_stream.build();
	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Building Output Stream 1");
	m_output1 = output_stream1.build(get_root(), "map1",  out[0]);
	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Building Output Stream 2");
	m_output2 = output_stream2.build(get_root(), "map2",  out[1]);

	// Do the actual work
	std::string key;
	std::string value;
	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Entering main loop");
	while (m_input->read(key, value)) {
		m_mapper->dual_map(key, value, *m_output1, *m_output2);
		m_num_proc++;
		send_update();
	}

	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Installing metadata");
	m_mapper->install_metadata1(out[0].metadata());
	m_mapper->install_metadata2(out[1].metadata());

	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Closing outputs");
	m_output1->close();
	m_output2->close();

	// SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Writing manifest");
	set_output(out);

	SPLOG_P(LOG_DEBUG, "dual_map_part_task::run> Done");
}

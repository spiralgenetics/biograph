#include "modules/mapred/splitter_task.h"
#include "modules/mapred/map_task.h"

REGISTER_TASK(splitter_task);


void splitter_task::prepare(manifest& input_manifest, double start)
{
	std::vector<input_stream_params> the_in_stream_params;
	m_split = manifest(input.get_sort(), 1);
	input_manifest.split_by_splitter(m_split, the_in_stream_params, splitter);

	SPLOG("splitter_task::prepare: %lu split subtasks needed, %lu files are already split correctly.",
		the_in_stream_params.size(), m_split.get_num_records());

	if (the_in_stream_params.size() == 0)
	{
		set_output(m_split);
		return;
	}

	for(unsigned int i = 0; i < the_in_stream_params.size(); i++)
	{
		std::unique_ptr<map_part_task> the_part_task = make_unique<map_part_task>();
		
		the_part_task->input_stream = the_in_stream_params[i];
		the_part_task->output_stream.num_partitions = num_partitions;
		the_part_task->output_stream.sort = input_manifest.get_sort();
		the_part_task->output_stream.split = splitter;
		the_part_task->output_stream.begin_on = the_in_stream_params[i].begin_on;
		the_part_task->map = "identity";
		the_part_task->map_param = map_param;
		the_part_task->update_freq = update_freq;

		m_subtasks.push_back(add_subtask(std::move(the_part_task)));
		update_progress(start + (1.0 - start)*i/input.get_size());
	}
}


void splitter_task::run()
{
	SPLOG("splitter_task::run started.");
	if (m_subtasks.size() == 0 && input.get_num_records() != 0)
	{
		prepare(input, 0.0);
		return;
	}

	manifest output_manifest(input.get_sort(), num_partitions);

	for(unsigned int i = 0; i < m_subtasks.size(); i++)
	{
		manifest subtask_output_manifest;
		get_output(subtask_output_manifest, m_subtasks[i]);
		output_manifest.add(subtask_output_manifest);
		update_progress(.85 * static_cast<double>(i) / static_cast<double>(m_subtasks.size()));
	}
	output_manifest.sort_file_infos();
	set_output(output_manifest);
}

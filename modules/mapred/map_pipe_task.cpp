#include <exception>

#include <boost/any.hpp>
#include <boost/filesystem.hpp>

#include "base/base.h"
#include "modules/mapred/ex_im_porter_data.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/bio_format/exporter.h"
#include "modules/bio_format/importer.h"
#include "modules/mapred/map_pipe_task.h"
#include "modules/mapred/pipe_mapper.h"
#include "modules/mapred/pipe_params.h"
#include "modules/mapred/unix_pipeline.h"

REGISTER_TASK(map_pipe_task);

subtask_id map_task::make_map_pipe_task(const input_stream_params & the_input_stream_params)
{
	auto task = make_unique<map_pipe_task>();

	task->input_stream = the_input_stream_params;
	task->output_stream.goal_size = output_goal_size;
	task->output_stream.num_partitions = num_partitions;
	if (stable_sort && input.get_sort() != "" && sort == "")
	{
		task->output_stream.presorted = true;
		task->output_stream.sort = input.get_sort();
	}
	else
	{
		task->output_stream.sort = sort;
	}

	task->output_stream.reduce = reduce;
	task->output_stream.reduce_param = reduce_param;
	task->map = map;
	task->map_param = map_param;
	task->update_freq = update_freq;

	return add_subtask(std::move(task));
}


void map_pipe_task::run()
{

	manifest output_manifest;
	m_input = input_stream.build();
	m_output = output_stream.build(get_root(), "map",  output_manifest);

	pipe_params the_pipe_params;
	json_deserialize(the_pipe_params, map_param);

	std::vector<std::string> temp_file_paths = create_temp_files(the_pipe_params);

	if (the_pipe_params.are_keys_sorted)
	{
		ex_im_porter_data the_ex_im_porter_data;
		msgpack_deserialize(the_ex_im_porter_data, the_pipe_params.ex_im_porter_data);
		the_ex_im_porter_data.start_key = input_stream.begin_on;
		the_ex_im_porter_data.end_key = input_stream.end_before;
		the_pipe_params.ex_im_porter_data = msgpack_serialize(the_ex_im_porter_data);
	}

	pipe_mapper_buffer the_mapper_buffer(*m_input, this);
	unix_pipeline a_pipe_writer(the_mapper_buffer, the_pipe_params.command, the_pipe_params.args, the_pipe_params.working_dir, [this]() {this->keep_alive();} );
	std::unique_ptr<importer> pipe_importer(
		importer_registry
		::get_safe(the_pipe_params.importer_type, the_mapper_buffer, true, the_pipe_params.ex_im_porter_data));
	std::unique_ptr<exporter> pipe_exporter(
		exporter_registry
		::get_safe(the_pipe_params.exporter_type, a_pipe_writer, true, the_pipe_params.ex_im_porter_data));

	try
	{
		pipe_mapper a_pipe_mapper(the_mapper_buffer, *pipe_exporter, *pipe_importer, a_pipe_writer);
		a_pipe_mapper(*m_output);
	}
	catch(...)
	{
		a_pipe_writer.log_child_stderr();
		std::for_each(temp_file_paths.begin(), temp_file_paths.end(),
			[](const std::string& the_path) { boost::filesystem::remove(the_path); }
		);
		throw;
	}

	m_output->close();
	set_output(output_manifest);
	std::for_each(temp_file_paths.begin(), temp_file_paths.end(),
		[](const std::string& the_path) { boost::filesystem::remove(the_path); }
	);

}


void map_pipe_task::processed_a_record()
{
	++m_records_processed;
	keep_alive();
}


void map_pipe_task::keep_alive() const
{
	update_progress(static_cast<double>(m_records_processed) / input_stream.num_records);
}


std::vector<std::string> map_pipe_task::create_temp_files(pipe_params& the_pipe_params) const
{
	std::vector<std::string> temp_file_paths;
	temp_file_paths.reserve(the_pipe_params.temp_files.size());

	for (auto temp_file_spec : the_pipe_params.temp_files)
	{
		std::string temp_file_path = create_temp_file(temp_file_spec);
		try
		{
			CHECK(the_pipe_params.args[temp_file_spec.arg_index].empty());
			the_pipe_params.args.at(temp_file_spec.arg_index) = temp_file_path;
		}
		catch(std::exception& e)
		{
			std::ostringstream error_stream;
			error_stream << "A temp file for command "
				<< the_pipe_params.command
				<< " was requested at argument "
				<< temp_file_spec.arg_index
				<< ", but there are only "
				<< the_pipe_params.args.size()
				<< " arguments: "
				<< e.what();
			throw io_exception(error_stream.str());
		}
		temp_file_paths.push_back(temp_file_path);
	}

	return temp_file_paths;
}


std::string map_pipe_task::create_temp_file(const temp_file_spec& the_temp_file_spec) const
{
	std::string temp_file_path = ::path(CONF_S(temp_root) + "/spiral-XXXXXX").bare_path();
	SPLOG("map_pipe_task::create_temp_file creating temp file %s", temp_file_path.c_str());

	file_writer temp_file_writer(temp_file_path);
	manifest_reader temp_manifest_reader(the_temp_file_spec.data);

	std::unique_ptr<exporter> temp_file_exporter_ptr(
		registry<exporter, writable&, bool, const std::string&>
		::get_safe(the_temp_file_spec.exporter_type, temp_file_writer, true, the_temp_file_spec.ex_im_porter_data));
	temp_file_writer.close();

	return temp_file_path;
}

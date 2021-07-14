#include <exception>

#include "modules/io/file_io.h"
#include "modules/mapred/run_wrapper_task.h"
#include "modules/mapred/pipe_params.h"
#include "modules/mapred/unix_pipeline.h"

REGISTER_TASK(run_wrapper_task);

void run_wrapper_task::run()
{
	null_writable devnull;

	path(params.working_dir).mkdir();

	for(size_t i = 0; i < inputs.size(); i++) {
		file_writer fw(params.working_dir + "/in_" + std::to_string(i));
		std::string json = json_serialize(inputs[i]);
		fw.write(json.c_str(), json.size());
		fw.close();
	}

	unix_pipeline the_pipe(devnull, params.command, params.args, params.working_dir, [this]() {this->keep_alive();} );

	the_pipe.close();

	std::vector<manifest> manifests(num_outputs);
	for(size_t i = 0; i < num_outputs; i++) {
		try {
			json_deserialize(manifests[i], slurp_file(params.working_dir + "/out_" + std::to_string(i)));
		}
		catch(const io_exception& e) {
			throw(io_exception(printstring("Task created invalid output out_%lu", i)));
		}
	}
	set_output(manifests);
}

void run_wrapper_task::keep_alive() const
{
	if(!update_progress(0.0)) {
		throw io_exception("Task is no longer in the run state, aborting.");
	};
}

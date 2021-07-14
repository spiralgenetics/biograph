#pragma once

#include "modules/mapred/map_task.h"
#include "modules/io/log.h"

struct pipe_params;
struct temp_file_spec;
class map_pipe_task : public task_impl<map_pipe_task>
{
public:
    static std::string s_type() { return "map_pipe"; }

    TRANSFER_OBJECT
    {
        VERSION(0);
        FIELD(input_stream, TF_STRICT);
        FIELD(output_stream, TF_STRICT);
        FIELD(map, TF_STRICT);
        FIELD(map_param, TF_STRICT);
        FIELD(update_freq, TF_STRICT);
    };

    void run();

    input_stream_params input_stream;
    output_stream_params output_stream;

    std::string map;
    std::string map_param;

    size_t update_freq;

	void processed_a_record();
	void keep_alive() const;

private:
    std::unique_ptr<kv_source> m_input;
    std::unique_ptr<kv_sink> m_output;

	unsigned int m_records_processed;

	// Creates temp files from the mainfests in the pipe params and sets up the command line arguments.
	std::vector<std::string> create_temp_files(pipe_params& the_pipe_params) const;
	std::string create_temp_file(const temp_file_spec& the_temp_file_spec) const;
};

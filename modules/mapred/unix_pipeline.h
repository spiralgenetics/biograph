#pragma once

#include "modules/io/io.h"

#include <vector>
#include <string>

#include <boost/circular_buffer.hpp>

#include <unistd.h>

class child_signalled : public io_exception
{
public:
	child_signalled(int wait_status)
		: io_exception("Child process was terminated by a signal")
		, m_signal(WTERMSIG(wait_status))
		{}
	int get_signal() const { return m_signal; }

private:
	int m_signal;
};

class map_pipe_task;
class unix_pipeline : public writable
{
public:
	unix_pipeline(
		// The object to which the data processed by the child will be written.
		// The object must stay alive until this object is destroyed.
		writable & processed_output_dest,
		// The path to the file to be executed.
		const std::string& command,
		// A vector of the arguments to the command.  Does not include the file to be executed.
		const std::vector<std::string>& arguments = std::vector<std::string>(),
		// Change the command current working directory.  If blank, the caller's value will be used.
		const std::string& working_dir = std::string(),
		// Optional function that will be updated (default: do nothing)
		const std::function<void()>& callback = std::function<void()>()
	);

	~unix_pipeline();

	void write(const char* buf, size_t len) override;

	void close() override;

	const boost::circular_buffer<char> & get_error_buffer() const { return m_child_error_buffer; }
	void clear_error_buffer() { m_child_error_buffer.clear(); }
	void log_child_stderr();
	pid_t get_child_process_id() const { return m_child_process_id; }

	// Throws a child_signalled exception if the child was signalled.
	// Returns nonsense if called before close().
	int get_child_exit_code() const;
	bool is_child_alive() const { return m_is_child_alive; }

private:
	int m_p2c_pipe_read_fd; // Parent to child pipe.
	int m_p2c_pipe_write_fd;

	int m_c2p_pipe_read_fd; // Child to parent pipe.
	int m_c2p_pipe_write_fd;

	int m_error_pipe_read_fd; // Pipe for reporting errors from the child.
	int m_error_pipe_write_fd;

	std::string m_child_command;
	pid_t m_child_process_id;

	bool m_is_child_alive;
	bool m_is_stdout_pipe_open;
	bool m_is_stderr_pipe_open;

	boost::circular_buffer<char> m_child_error_buffer;
	static const size_t	mk_child_error_buffer_size = 16 * 1024;

	writable & m_processed_output_destination;

	int m_child_wait_status;

	std::function<void()> m_callback;

	void create_pipes();
	// The close_file_descriptor method will add errno and the error to the error message.
	void close_file_descriptor(int file_descriptor, const std::string& error_message);
	void move_file_descriptor(int source_file_descriptor, int dest_file_dscriptor, const std::string& error_message);
	void start_child(const std::string& command,
		const std::vector<std::string>& arguments,
		const std::string& working_dir
	);
	// Caution! The returned char pointers in the vector will be invalid when "command" and "argument"
	// go out of scope.
	std::vector<char*> build_argument_array(const std::string& command, const std::vector<std::string>& arguments);

	void read_child_error_stream(int error_file_descriptor);
	bool read_child_stdout(int stdout_file_descriptor); // Returns whether pipe is is still open
	size_t write_to_child(int child_write_file_descriptor, const char* buf, size_t len);
	void wait_for_child(bool is_blocking) throw();

	static const bool k_block = true;
	static const bool k_no_block = false;
};

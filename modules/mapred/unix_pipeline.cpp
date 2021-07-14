#include "modules/mapred/unix_pipeline.h"
#include "modules/io/log.h"
#include "modules/io/config.h"
#include "modules/mapred/map_pipe_task.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <cstring>
#include <algorithm>

#include <unistd.h>
#include <errno.h>
#include <sys/wait.h>
#include <fcntl.h>

unix_pipeline::unix_pipeline(
	writable & processed_output_dest,
	const std::string& command,
	const std::vector<std::string>& arguments,
	const std::string& working_dir_path,
	const std::function<void()>& callback
)
		: m_child_command(command)
		, m_is_child_alive(false)
		, m_is_stdout_pipe_open(false)
		, m_is_stderr_pipe_open(false)
		, m_child_error_buffer(mk_child_error_buffer_size)
		, m_processed_output_destination(processed_output_dest)
		, m_child_wait_status(0)
		, m_callback(callback)
{
	create_pipes();
	start_child(command, arguments, working_dir_path);
}

void unix_pipeline::write(const char* buf, size_t len)
{
	int max_file_descriptor = 0;

	fd_set read_file_descriptor_set;
	fd_set write_file_descriptor_set;

	while (len > 0) {
		FD_ZERO(&read_file_descriptor_set);
		if (m_is_stdout_pipe_open) {
			FD_SET(m_c2p_pipe_read_fd, &read_file_descriptor_set);
			max_file_descriptor = std::max(max_file_descriptor, m_c2p_pipe_read_fd);
		}
		if (m_is_stderr_pipe_open) {
			FD_SET(m_error_pipe_read_fd, &read_file_descriptor_set);
			max_file_descriptor = std::max(max_file_descriptor, m_error_pipe_read_fd);
		}

		FD_ZERO(&write_file_descriptor_set);
		FD_SET(m_p2c_pipe_write_fd, &write_file_descriptor_set);
		max_file_descriptor = std::max(max_file_descriptor, m_p2c_pipe_write_fd);

		struct timeval select_timeout = {30, 0};

		int select_ret = ::select(max_file_descriptor + 1,
			&read_file_descriptor_set,
			&write_file_descriptor_set,
			NULL, &select_timeout);

		// Update progress
		if (m_callback) {
			m_callback();
		}

		// Error
		if (select_ret == -1) {
			throw io_exception("unix_pipeline::write> Call to select failed");
		}
		// Timeout
		else if (select_ret == 0) {
			SPLOG("unix_pipeline::write> timeout!");
			continue;
		}

		// Anything on child stderr?
		if (FD_ISSET(m_error_pipe_read_fd, &read_file_descriptor_set)) {
			read_child_error_stream(m_error_pipe_read_fd);
			log_child_stderr();
		}
		// Anything on child stdout?
		if (FD_ISSET(m_c2p_pipe_read_fd, &read_file_descriptor_set)) {
			if (! read_child_stdout(m_c2p_pipe_read_fd)) {
				throw io_exception("unix_pipeline::write> Child stdout pipe closed prematurely.");
			}
		}
		// Can we write to the child?
		if (FD_ISSET(m_p2c_pipe_write_fd, &write_file_descriptor_set)) {
			ssize_t write_count = ::write(m_p2c_pipe_write_fd, buf, len);
			if (write_count == -1) {
				throw io_exception("unix_pipeline::write_to_child> Write to child failed");
			}
			buf += write_count;
			len -= write_count;
		}
	}
}

// The parent blocks after calling close()
void unix_pipeline::close()
{
	close_file_descriptor(m_p2c_pipe_write_fd, "unix_pipeline::close> Failed to close parent to child write on overall close.");
	// SPLOG("unix_pipeline::close> Parent closed write pipeline to child %d successfully.", m_child_process_id);

	int	max_file_descriptor = std::max(m_error_pipe_read_fd, m_c2p_pipe_read_fd) + 1;

	fd_set read_file_descriptor_set;

	while (m_is_stdout_pipe_open || m_is_stderr_pipe_open) {

		FD_ZERO(&read_file_descriptor_set);
		if (m_is_stdout_pipe_open) {
			FD_SET(m_c2p_pipe_read_fd, &read_file_descriptor_set);
		}
		if (m_is_stderr_pipe_open) {
			FD_SET(m_error_pipe_read_fd, &read_file_descriptor_set);
		}

		struct timeval select_timeout = {30, 0};
		int select_ret = ::select(max_file_descriptor,
			&read_file_descriptor_set,
			NULL, NULL, &select_timeout);

		// Update progress
		if (m_callback) {
			m_callback();
		}

		// Error
		if (select_ret == -1) {
			throw io_exception("unix_pipeline::close> Call to select on close failed");
		}
		// Time out
		else if (select_ret == 0) {
			continue;
		}

		// Anything on child stderr?
		if (FD_ISSET(m_error_pipe_read_fd, &read_file_descriptor_set)) {
			read_child_error_stream(m_error_pipe_read_fd);
		}
		// Anything on child stdout?
		if (FD_ISSET(m_c2p_pipe_read_fd, &read_file_descriptor_set)) {
			read_child_stdout(m_c2p_pipe_read_fd);
		}
	}

	m_processed_output_destination.close();
	wait_for_child(k_block);
}

void unix_pipeline::create_pipes()
{
	int pipe_file_descriptors[2];
	int pipe_ret = ::pipe(pipe_file_descriptors);
	if (pipe_ret != 0) {
		throw io_exception("unix_pipeline::create_pipes> Parent to child pipe failed to open");
	}

	m_p2c_pipe_read_fd = pipe_file_descriptors[0];
	m_p2c_pipe_write_fd = pipe_file_descriptors[1];

	pipe_ret = ::pipe(pipe_file_descriptors);
	if (pipe_ret != 0) {
		throw io_exception("unix_pipeline::create_pipes> Child to parent pipe failed to open");
	}

	m_c2p_pipe_read_fd = pipe_file_descriptors[0];
	m_c2p_pipe_write_fd = pipe_file_descriptors[1];

	pipe_ret = ::pipe(pipe_file_descriptors);
	if (pipe_ret != 0) {
		throw io_exception("unix_pipeline::create_pipes> Error pipe failed to open");
	}

	m_error_pipe_read_fd = pipe_file_descriptors[0];
	m_error_pipe_write_fd = pipe_file_descriptors[1];
	m_is_stdout_pipe_open = true;
	m_is_stderr_pipe_open = true;
}

void unix_pipeline::start_child(const std::string& command,
		const std::vector<std::string>& arguments,
		const std::string& working_dir_path)
{
	std::string command_line = command + std::string(" ") + boost::join(arguments, " ");
	SPLOG("unix_pipeline::start_child> About to run exec with command line '%s'", command_line.c_str());

	m_child_process_id = ::fork();

	// Fork fail
	if (m_child_process_id == -1) {
		throw io_exception("unix_pipeline::start_child> Fork failed");
	}
	// Child process
	else if (m_child_process_id == 0) {
		close_file_descriptor(m_p2c_pipe_write_fd, "unix_pipeline::start_child> Parent to child write pipe failed to close");
		close_file_descriptor(m_c2p_pipe_read_fd, "unix_pipeline::start_child> Child to parent read pipe failed to close");
		close_file_descriptor(m_error_pipe_read_fd, "unix_pipeline::start_child> Child to parent error pipe failed to close");
		move_file_descriptor(m_p2c_pipe_read_fd, STDIN_FILENO, "unix_pipeline::start_child> Connecting parent to child read end pipe to child stdin failed");
		move_file_descriptor(m_c2p_pipe_write_fd, STDOUT_FILENO, "unix_pipeline::start_child> Connecting child to parent write end pipe to child stdout failed");
		move_file_descriptor(m_error_pipe_write_fd, STDERR_FILENO, "unix_pipeline::start_child> Connecting child error pipe to child stderr failed");

		if (! working_dir_path.empty()) {
			int chdir_ret = ::chdir(working_dir_path.c_str());
			if (chdir_ret != 0) {
				throw io_exception("unix_pipeline::start_child> chdir failed");
			}
		}

		std::vector<char*> argument_array = build_argument_array(command, arguments);
		int exec_ret = ::execv(command.c_str(), &argument_array[0]);
		if (exec_ret != 0) {
			std::string error_string = "unix_pipeline::start_child> Call to execv failed. Command = ";
			error_string += command;
			std::vector<std::string>::const_iterator argument_iterator = arguments.begin();
			while (argument_iterator != arguments.end()) {
				error_string += " ";
				error_string += *argument_iterator++;
			}
			error_string += ", working dir path = ";
			error_string += working_dir_path;
			io_exception e(error_string);
			SPLOG("%s", e.message().c_str());
			::exit(errno);
		}
	}
	// Parent process
	else {
		m_is_child_alive = true;
		close_file_descriptor(m_p2c_pipe_read_fd, "unix_pipeline::start_child> Parent to child read pipe failed to close");
		close_file_descriptor(m_c2p_pipe_write_fd, "unix_pipeline::start_child> Child to parent write pipe failed to close");
		close_file_descriptor(m_error_pipe_write_fd, "unix_pipeline::start_child> Parent to child error pipe failed to close");

		// We must set close on exec for the parent pipes so that a future child in another
		// instance of this class doesn't inherit these pipes.
		fcntl(m_p2c_pipe_write_fd, F_SETFD, fcntl(m_p2c_pipe_write_fd, F_GETFD) | FD_CLOEXEC);
		fcntl(m_c2p_pipe_read_fd, F_SETFD, fcntl(m_c2p_pipe_read_fd, F_GETFD) | FD_CLOEXEC);
		fcntl(m_error_pipe_read_fd, F_SETFD, fcntl(m_error_pipe_read_fd, F_GETFD) | FD_CLOEXEC);
	}
}

std::vector<char*> unix_pipeline::build_argument_array(const std::string& command, const std::vector<std::string>& arguments)
{
	std::vector<char*> argument_array;
	argument_array.reserve(arguments.size() + 2);
	argument_array.push_back(const_cast<char*>(command.c_str()));

	std::vector<std::string>::const_iterator argument_iterator = arguments.begin();
	while (argument_iterator != arguments.end()) {
		argument_array.push_back(const_cast<char*>(argument_iterator->c_str()));
		argument_iterator++;
	}
	argument_array.push_back(NULL);

	return argument_array;
}

void unix_pipeline::close_file_descriptor(int file_descriptor, const std::string& error_message)
{
	int close_ret = ::close(file_descriptor);
	if (close_ret != 0) {
		throw io_exception(error_message);
	}
}

void unix_pipeline::move_file_descriptor(int source_file_descriptor, int dest_file_descriptor, const std::string& error_message)
{
	int dup2_ret = ::dup2(source_file_descriptor, dest_file_descriptor);
	if (dup2_ret == -1) {
		throw io_exception(error_message);
	}
	std::string close_error = std::string("unix_pipeline::move_file_descriptor> close error: ") + error_message;
	close_file_descriptor(source_file_descriptor, error_message);
}

void unix_pipeline::read_child_error_stream(int error_file_descriptor)
{
	if (! m_is_stderr_pipe_open) {
		return;
	}

	std::vector<char> read_buffer(PIPE_BUF);
	ssize_t	read_count = ::read(error_file_descriptor, &read_buffer[0], read_buffer.size());
	read_buffer.resize(read_count);
	// Error
	if (read_count == -1) {
		throw io_exception("unix_pipeline::read_child_error_stream> Child error stream read failed");
	}
	// End of file, i.e. child closed stdout.
	else if (read_count == 0) {
		close_file_descriptor(error_file_descriptor, "unix_pipeline::read_child_error_stream> Child error pipe failed to close");
		m_is_stderr_pipe_open = false;
		// Check if child crashed.
		wait_for_child(k_no_block);
	}
	// We got something
	else {
		std::copy(read_buffer.begin(), read_buffer.begin() + read_count, std::back_inserter(m_child_error_buffer));
	}
}

bool unix_pipeline::read_child_stdout(int stdout_file_descriptor)
{
	if (! m_is_stdout_pipe_open) {
		return false;
	}

	std::vector<char> read_buffer(PIPE_BUF);
	ssize_t	read_count = ::read(stdout_file_descriptor, &read_buffer[0], read_buffer.size());
	read_buffer.resize(read_count);

	// Error
	if (read_count == -1) {
		throw io_exception("unix_pipeline::read_child_stdout> Child output stream read failed");
	}
	// End of file, i.e. child closed stdout.
	else if (read_count == 0) {
		m_is_stdout_pipe_open = false;
		close_file_descriptor(stdout_file_descriptor, "unix_pipeline::read_child_stdout> Child stdout pipe failed to close");
		// Check if child crashed
		wait_for_child(k_no_block);
	}
	// We got some processed output, pass it on...
	else {
		m_processed_output_destination.write(&read_buffer[0], read_count);
	}

	return m_is_stdout_pipe_open;
}

void unix_pipeline::wait_for_child(bool is_blocking) throw()
{
	if (! m_is_child_alive) {
		return;
	}

	pid_t wait_ret = ::waitpid(m_child_process_id, &m_child_wait_status, WNOHANG);
	if (wait_ret == 0 && is_blocking) {
		::sleep(1);
		wait_ret = ::waitpid(m_child_process_id, &m_child_wait_status, 0);
	}

	if (wait_ret == -1) {
		io_exception e("unix_pipeline::wait_for_child> Waiting for child failed");
		SPLOG("%s", e.message().c_str());
	}
	else if (wait_ret == m_child_process_id) {
		if (WIFEXITED(m_child_wait_status)) {
			SPLOG("unix_pipeline::wait_for_child> Child process %s with pid %d exited with status %d",
				m_child_command.c_str(), m_child_process_id, WEXITSTATUS(m_child_wait_status));
				m_is_child_alive = false;
		}
		else if (WIFSIGNALED(m_child_wait_status)) {
			SPLOG("unix_pipeline::wait_for_child> Child process %s with pid %d was signalled with signal %d",
				m_child_command.c_str(), m_child_process_id, WTERMSIG(m_child_wait_status));
			m_is_child_alive = false;
		}
		else {
			SPLOG("unix_pipeline::wait_for_child> returned an unknown status %d", m_child_wait_status);
		}
	}
	else if (! is_blocking && wait_ret == 0) {
		// Do nothing.
	}
	else {
		std::string error_string = "Unexpected return from waiting for child of ";
		error_string += boost::lexical_cast<std::string>(wait_ret);
		SPLOG("unix_pipeline::wait_for_child> %s", error_string.c_str());
	}
}

int unix_pipeline::get_child_exit_code() const
{
	if (m_is_child_alive) {
		throw io_exception("unix_pipeline::get_child_exit_code> Cannot get child exit code until child has exited!");
	}
	else if (WIFSIGNALED(m_child_wait_status)) {
		throw child_signalled(m_child_wait_status);
	}
	else if (! WIFEXITED(m_child_wait_status)) {
		throw io_exception("unix_pipeline::get_child_exit_code> Child did not exit normally, but no signal was detected.");
	}

	return WEXITSTATUS(m_child_wait_status);
}

unix_pipeline::~unix_pipeline()
{
	if (m_is_child_alive) {
		wait_for_child(k_no_block);
		::sleep(1);
		wait_for_child(k_no_block);
		int kill_ret = ::kill(m_child_process_id, SIGTERM);
		SPLOG("~unix_pipeline> Sending child %d SIGTERM.", m_child_process_id);
		if (kill_ret == -1) {
			io_exception e("Attempt to terminate child process failed");
			SPLOG("~unix_pipeline> %s", e.message().c_str());
		}
		::sleep(3);
		wait_for_child(k_no_block);
		// There's a bit of a race condition here if the child dies between the
		// kill 0 and the kill SIGKILL, but nothing beyond logging an extraneous
		// error willl occur.
		if (::kill(m_child_process_id, 0) == 0) {
			SPLOG("~unix_pipeline> Sending child %d SIGKILL.", m_child_process_id);
			kill_ret = ::kill(m_child_process_id, SIGKILL);
			if (kill_ret == -1) {
				io_exception e("Attempt to kill child process failed");
				SPLOG("~unix_pipeline> %s", e.message().c_str());
			}
			::sleep(1);
			wait_for_child(k_no_block);
		}
	}
}

void unix_pipeline::log_child_stderr()
{
	SPLOG("unix_pipeline::log_child_stderr> Logging stderr from child %d", m_child_process_id);
	auto child_stderr = get_error_buffer();
	if (! child_stderr.empty()) {
		// This is wasteful but necessary because the circular buffer iterators are broken.
		std::string stderr_copy;
		stderr_copy.reserve(child_stderr.size());
		std::copy(child_stderr.begin(), child_stderr.end(), std::back_inserter(stderr_copy));
		boost::char_separator<char> newline_separator("\n");
		boost::tokenizer<boost::char_separator<char>> the_lines(stderr_copy, newline_separator);
		auto line_iter = the_lines.begin();
		while (line_iter != the_lines.end()) {
			SPLOG("unix_pipeline::log_child_stderr> %s", line_iter++->c_str());
		}
	}
}

#include "modules/io/command.h"
#include "modules/io/utils.h"
#include "modules/io/make_unique.h"
#include "modules/io/log.h"
#include "modules/io/mem_io.h"


#include <boost/filesystem.hpp>

#include <future>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <unistd.h>
#include <sys/wait.h>

namespace exec
{

pipe_reader::pipe_reader(int fd)
	: m_fd(fd)
{}

pipe_reader::~pipe_reader()
{
	close();
}

void pipe_reader::close()
{
	// SPLOG("pipe_reader::close>");
	if (m_fd != -1) {
		if (::close(m_fd)) {
			throw io_exception(printstring("exec::pipe_reader::close> ::close() failed: %s",
				strerror(errno)));
		}
		m_fd = -1;
	}
}

int pipe_reader::fileno() const
{
	return m_fd;
}

int pipe_reader::base_read(char* buf, size_t len)
{
	return ::read(m_fd, buf, len);
}

pipe_writer::pipe_writer(int fd)
	: m_fd(fd)
{}

pipe_writer::~pipe_writer()
{
	close();
}

int pipe_writer::fileno() const
{
	return m_fd;
}

int pipe_writer::base_write(const char* buf, int len)
{
	return ::write(m_fd, buf, len);
}

int pipe_writer::base_flush()
{
	return ::fsync(m_fd);
}

int pipe_writer::base_close()
{
	// SPLOG("pipe_writer::base_close>");
	int retval = 0;
	if (m_fd != -1) {
		retval = ::close(m_fd);
		if (!retval) {
			m_fd = -1;
		}
	}
	return retval;
}

pipe::pipe()
{
	int fds[] = { -1, -1 };
	auto retval = ::pipe(fds);
	if (retval) {
		throw io_exception(printstring("exec::pipe> ::pipe() failed: %s",
			strerror(errno)));
	}
	m_reader = make_unique<pipe_reader>(fds[0]);
	m_writer = make_unique<pipe_writer>(fds[1]);
}

void pipe::close_read()
{
	m_reader.reset();
}

void pipe::close_write()
{
	m_writer.reset();
}

void pipe::dup_read(int target)
{
	if (!m_reader) {
		throw io_exception("exec::pipe::dup_read> Invalid state, m_reader is nil");
	}
	if (::dup2(m_reader->fileno(), target) == -1) {
		throw io_exception(printstring("exec::pipe::dup_read> ::dup2() failed: %s",
			strerror(errno)));
	}
}

void pipe::dup_write(int target)
{
	if (!m_writer) {
		throw io_exception("exec::pipe::dup_write> Invalid state, m_writer is nil");
	}
	if (::dup2(m_writer->fileno(), target) == -1) {
		throw io_exception(printstring("exec::pipe::dup_write> ::dup2() failed: %s",
			strerror(errno)));
	}
}

readable* pipe::reader()
{
	return m_reader.get();
}

writable* pipe::writer()
{
	return m_writer.get();
}

command::command(const std::string& path, const std::vector<std::string>& args)
	: m_path(path)
	, m_args(args)
{}

int command::run()
{
	start();
	return wait();
}

void command::start()
{
	std::string command_line_string{m_path + " "};
	for (const auto& arg : m_args)
	{
		command_line_string += '\'';
		command_line_string += arg;
		command_line_string += '\'';
		command_line_string += ' ';
	}
	SPLOG("command::start> %s", command_line_string.c_str());

	m_pid = fork();
	switch (m_pid) {
	case -1: // failed to fork
		on_error();
		break;
	case 0: // child
		on_child();
		break;
	default: // parent
		on_parent();
		break;
	}
}

int command::wait()
{
	int status = 0;
	while (true) {
		int retcode = waitpid(m_pid, &status, 0);
		if (retcode != -1 || errno != EINTR) {
			break;
		}
	}

	m_pid = -1;
	if (WIFSIGNALED(status)) {
		return WTERMSIG(status);
	}
	else if (WIFEXITED(status)) {
		return WEXITSTATUS(status);
	}
	return 0;
}

writable* command::stdin()
{
	if (!m_stdin) {
		m_stdin = make_unique<pipe>();
	}
	return m_stdin->writer();
}

readable* command::stdout()
{
	if (!m_stdout) {
		m_stdout = make_unique<pipe>();
	}
	return m_stdout->reader();
}

readable* command::stderr()
{
	if (!m_stderr) {
		m_stderr = make_unique<pipe>();
	}
	return m_stderr->reader();
}

std::string command::path() const
{
	return m_path;
}

std::vector<std::string> command::args() const
{
	return m_args;
}

void command::on_error()
{
	m_stdin.reset();
	m_stdout.reset();
	m_stderr.reset();
	throw io_exception(printstring("exec::command::on_error> fork() failed: %s",
		strerror(errno)));
}

void command::on_child()
{
	try {
		// NOTE: use only 'async-signal-safe' system functions between here and the exec*()
		// http://pubs.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_04.html#tag_02_04_03
		if (m_stdin) {
			m_stdin->close_write();
			m_stdin->dup_read(STDIN_FILENO);
		}
		if (m_stdout) {
			m_stdout->close_read();
			m_stdout->dup_write(STDOUT_FILENO);
		}
		if (m_stderr) {
			m_stderr->close_read();
			m_stderr->dup_write(STDERR_FILENO);
		}

		std::vector<const char*> args;
		args.push_back(m_path.c_str());
		for (const auto& arg : m_args) {
			args.push_back(arg.c_str());
		}
		args.push_back(nullptr);

		auto const_args = const_cast<char* const*>(args.data());
		(void) execvp(m_path.c_str(), const_args);
		// only way we fall thru to here is if execv() failed
		std::exit(errno);
	}
	catch (const io_exception& io) {
		// NOTE: unfortunately we can't use syslog in this context
		// see: http://pubs.opengroup.org/onlinepubs/009695399/functions/pthread_atfork.html
		// also: https://sourceware.org/ml/libc-alpha/2006-09/msg00045.html
		fprintf(::stderr, "Uncaught exception in forked child: %s\n", io.message().c_str());
		std::exit(1);
	}
}

void command::on_parent()
{
	if (m_stdin) {
		m_stdin->close_read();
	}
	if (m_stdout) {
		m_stdout->close_write();
	}
	if (m_stderr) {
		m_stderr->close_write();
	}
}

int call(const std::string& path, const std::vector<std::string>& args)
{
	exec::command cmd(path, args);
	return cmd.run();
}

void check_call(const std::string& path, const std::vector<std::string>& args)
{
  mem_io errbuf("", track_alloc("command:check_call"));
	exec::command cmd(path, args);
	auto stderr = cmd.stderr();
	cmd.start();
	io_copy(*stderr, errbuf);
	auto retcode = cmd.wait();
	if (retcode) {
		auto msg = printstring("check_call> %s failed with retcode: %d", path.c_str(), retcode);
		SPLOG("%s", msg.c_str());
		SPLOG("%s", errbuf.str().c_str());
		throw io_exception(msg);
	}
}

std::string check_output(const std::string& path, const std::vector<std::string>& args)
{
  mem_io outbuf("",track_alloc("check_output:outbuf"));
  mem_io errbuf("",track_alloc("check_outpout:errbuf"));
	exec::command cmd(path, args);
	auto stdout = cmd.stdout();
	auto stderr = cmd.stderr();

	cmd.start();

	// copy stdout and stderr simultanously so that neither run out of buffer space
	// if either pipe runs out of memory, io_copy() will hang forever
	io_copy_pairs({
		{stdout, &outbuf},
		{stderr, &errbuf}
	});

	auto retcode = cmd.wait();
	if (retcode) {
		auto msg = printstring("check_output> %s failed with retcode: %d", path.c_str(), retcode);
		SPLOG("%s", msg.c_str());
		SPLOG("%s", errbuf.str().c_str());
		throw io_exception(msg);
	}

	// SPLOG("%s", errbuf.str().c_str());
	return outbuf.str();
}

std::string communicate(
	readable& reader,
	const std::string& path,
	const std::vector<std::string>& args)
{
  mem_io outbuf("",track_alloc("communicate:outbuf"));
  mem_io errbuf("",track_alloc("communicate:errbuf"));
	exec::command cmd(path, args);
	auto stdin = cmd.stdin();
	auto stdout = cmd.stdout();
	auto stderr = cmd.stderr();

	cmd.start();

	io_copy(reader, *stdin);
	stdin->close();

	// copy stdout and stderr simultanously so that neither run out of buffer space
	// if either pipe runs out of memory, io_copy() will hang forever
	io_copy_pairs({
		{stdout, &outbuf},
		{stderr, &errbuf}
	});

	auto retcode = cmd.wait();
	if (retcode) {
		auto msg = printstring("communicate> %s failed with retcode: %d", path.c_str(), retcode);
		SPLOG("%s", msg.c_str());
		SPLOG("%s", errbuf.str().c_str());
		throw io_exception(msg);
	}

	// SPLOG("%s", errbuf.str().c_str());
	return outbuf.str();
}

int ignore_io(
	const std::string& path,
	const std::vector<std::string>& args)
{
	null_writable outbuf;
	null_writable errbuf;
	exec::command cmd(path, args);
	auto stdin = cmd.stdin();
	auto stdout = cmd.stdout();
	auto stderr = cmd.stderr();

	cmd.start();

	stdin->close();

	io_copy_pairs({
		{stdout, &outbuf},
		{stderr, &errbuf}
	});

	return cmd.wait();
}

std::string get_exe_dir(pid_t process_id)
{
	std::string kernel_proc_path = "/proc/" + std::to_string(process_id) + "/exe";
	std::array<char, 4096> exe_path;
	ssize_t exe_path_len = ::readlink(kernel_proc_path.c_str(), exe_path.data(), exe_path.size());

	if (exe_path_len == -1) {
		int saved_errno = errno;
		std::array<char, 1024> error_string;
		char* error_buf = ::strerror_r(saved_errno, error_string.data(), error_string.size());
		if (error_buf) {
			throw io_exception(boost::format("exec::get_exe_dir> readlink failed with proc path \"%1%\", process ID %2% and errno %3%, \"%4%\"")
				% kernel_proc_path % process_id % saved_errno % error_string.data());
		} else {
			throw io_exception(boost::format("exec::get_exe_dir> strerror failed on errno = %1%") % saved_errno);
		}
	} else if (exe_path_len > static_cast<ssize_t>(exe_path.size() - 1)) {
		throw io_exception(boost::format("exec::get_exe_dir> readlink buffer is too small.  Had %1% bytes, but tried to write %2% bytes.")
			% exe_path.size() % exe_path_len);
	}

	exe_path[exe_path_len] = '\0';

	return boost::filesystem::path(exe_path.data()).parent_path().native();
}

} // namespace exec

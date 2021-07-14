#pragma once

#include "modules/io/io.h"
#include <memory>
#include <vector>
#include <string>

namespace exec
{

class pipe_reader : public read_wrapper
{
public:
	pipe_reader(int fd);
	~pipe_reader();

	void close() override;
	int fileno() const;

private:
	int base_read(char* buf, size_t len) override;

private:
	int m_fd = -1;
};

class pipe_writer : public write_wrapper
{
public:
	pipe_writer(int fd);
	~pipe_writer();

	int fileno() const;

private:
	int base_write(const char* buf, int len) override;
	int base_flush() override;
	int base_close() override;

	int m_fd = -1;
};

class pipe
{
public:
	pipe();

	void close_read();
	void close_write();

	void dup_read(int target);
	void dup_write(int target);

	readable* reader();
	writable* writer();

private:
	std::unique_ptr<pipe_reader> m_reader;
	std::unique_ptr<pipe_writer> m_writer;
};

class command
{
public:
	command(const command&) = delete;
	command& operator =(const command&) = delete;

	// Do not pass the path as the zeroeth argument, the class will do that later.
	command(const std::string& path, const std::vector<std::string>& args = {});

	int run();
	void start();
	int wait();

	writable* stdin();
	readable* stdout();
	readable* stderr();

	std::string path() const;
	std::vector<std::string> args() const;

private:
	void on_error();
	void on_child();
	void on_parent();

private:
	std::string m_path;
	std::vector<std::string> m_args;
	std::unique_ptr<pipe> m_stdin;
	std::unique_ptr<pipe> m_stdout;
	std::unique_ptr<pipe> m_stderr;
	pid_t m_pid = -1;
};

// Do not pass the path as the zeroeth argument
int call(const std::string& path, const std::vector<std::string>& args = {});
void check_call(const std::string& path, const std::vector<std::string>& args = {});
std::string check_output(const std::string& path, const std::vector<std::string>& args = {});

// Pass a readable into stdin, log stderr and return stdout in a string.
// Throw an exception if a non-zero return code is encountered.
std::string communicate(
	readable& reader,
	const std::string& path, 
	const std::vector<std::string>& args = {});

// Close stdin and throw stdout and stderr away.  Return the return code.
int ignore_io(
	const std::string& path,
	const std::vector<std::string>& args);

// Get the executable directory of the given process
std::string get_exe_dir(pid_t process_id);

} // namespace exec

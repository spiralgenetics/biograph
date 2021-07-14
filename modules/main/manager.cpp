#include "modules/main/main.h"

#include "modules/io/log.h"
#include "modules/io/io.h"
#include "modules/io/utils.h"
#include "modules/io/config.h"
#include "base/command_line.h"

#include <thread>
#include <fcntl.h>
#include <sys/wait.h>
#include <sys/stat.h>

enum PipeSide
{
	READ,
	WRITE,
};

enum WorkerState
{
	RUNNING,
	OK,
	FAIL,
	NO_WORK,
	UPDATE,
};

struct Pipe
{
	Pipe()
		: fds{{ -1, -1 }}
	{
		if (::pipe(fds.data())) {
			throw io_exception(printstring("Pipe> ::pipe() failed: %s", strerror(errno)));
		}
	}

	~Pipe()
	{
		close(PipeSide::READ);
		close(PipeSide::WRITE);
	}

	void close(PipeSide side)
	{
		if (fds[side] != -1) {
			::close(fds[side]);
			fds[side] = -1;
		}
	}

	int side(PipeSide side)
	{
		return fds[side];
	}

	std::array<int, 2> fds;
};

class WorkerProcess
{
public:
	WorkerProcess(const std::string& profile);
	WorkerState poll();

private:
	int wait(bool hang);
	void terminate();

private:
	std::string m_profile;
	pid_t m_pid;
	Pipe m_pipe;
	std::chrono::steady_clock::time_point m_last_heard;
	std::chrono::seconds WORKER_TIMEOUT{CONF_T(int, task_timeout)};
};

class ManagerMain : public Main
{
protected:
	void add_args() override;
	int run(po::variables_map vars) override;

private:
	size_t m_num_procs = 0;
	std::string m_profile;
};

std::unique_ptr<Main> manager_main()
{
	return std::unique_ptr<Main>(new ManagerMain);
}

void touch(const std::string& path)
{
	int fd = ::open(path.c_str(), O_WRONLY|O_CREAT|O_NOCTTY|O_NONBLOCK, 0666);
	if (fd < 0) {
		SPLOG("touch> ::open() failed: %s", strerror(errno));
		return;
	}

	::close(fd);

	int retcode = ::utimensat(AT_FDCWD, path.c_str(), nullptr, 0);
	if (retcode) {
		SPLOG("touch> ::utimensat() failed: %s", strerror(errno));
	}
}

void ManagerMain::add_args()
{
	auto cpus = std::thread::hardware_concurrency();
	auto cpus_env = getenv("SPIRAL_NUM_WORKERS");
	if (cpus_env) {
		try {
			cpus = std::stoi(cpus_env);
		}
		catch (const std::exception& ex) {
			SPLOG("Ignoring invalid value for environment variable: SPIRAL_NUM_WORKERS.");
		}
	}

	std::string profile;
	if (auto c_profile = getenv("WORKER_PROFILE")) {
		profile = c_profile;
	}

	m_options.add_options()
		("num_procs", po::value(&m_num_procs)->default_value(cpus))
		("profile", po::value(&m_profile)->default_value(profile))
	;
}

int ManagerMain::run(po::variables_map vars)
{
	SPLOG("Starting manager: profile '%s', num_procs %zu", m_profile.c_str(), m_num_procs);

	typedef std::unique_ptr<WorkerProcess> WorkerProcess_ptr;
	std::vector<WorkerProcess_ptr> workers(m_num_procs);

	while (true) {
		for (auto& worker : workers) {
			if (worker) {
				auto state = worker->poll();
				if (state != WorkerState::UPDATE && state != WorkerState::RUNNING) {
					worker = nullptr;
				}
			}
		}

		// find the first empty slot
		auto it = std::find(workers.begin(), workers.end(), nullptr);
		if (it != workers.end()) {
			try {
				*it = WorkerProcess_ptr(new WorkerProcess(m_profile));
			}
			catch (const io_exception& ex) {
				SPLOG("Error: %s", ex.message().c_str());
			}
		}

		std::this_thread::sleep_for(std::chrono::milliseconds(1500));
	}

	return 0;
}

WorkerProcess::WorkerProcess(const std::string& profile)
	: m_profile(profile)
{
	m_pid = fork();
	switch (m_pid) {
	case -1: // failed to fork
		throw io_exception(printstring("WorkerProcess> ::fork() failed: %s", strerror(errno)));
		break;
	case 0: // child
		setproctitle("biograph_worker");
		log_change_name("biograph_worker");
		m_pipe.close(PipeSide::READ);
		::dup2(m_pipe.side(PipeSide::WRITE), STDOUT_FILENO);
		try {
			std::exit(do_worker(m_profile));
		}
		catch (const io_exception& ex) {
			SPLOG("Error: %s", ex.message().c_str());
			std::exit(1);
		}
		catch (...) {
			std::exit(1);
		}
		break;
	default: // parent
		m_pipe.close(PipeSide::WRITE);
		if (::fcntl(m_pipe.side(PipeSide::READ), F_SETFL, O_NONBLOCK) == -1) {
			throw io_exception(printstring("WorkerProcess> ::fcntl() failed: %s", strerror(errno)));
		}
		break;
	}

	m_last_heard = std::chrono::steady_clock::now();
}

WorkerState WorkerProcess::poll()
{
	char last_char = '\0';
	while (true) {
		auto result = ::read(m_pipe.side(PipeSide::READ), &last_char, 1);
		if (result == -1) {
			if (errno == EAGAIN) {
				break;
			}
			SPLOG("Could not read from process, pid: %d", m_pid);
			terminate();
			return WorkerState::FAIL;
		}
		if (result == 0) {
			break;
		}
		m_last_heard = std::chrono::steady_clock::now();
	}

	// determine disposition of completed worker
	if (wait(false)) {
		switch (last_char) {
		case 'S':
			return WorkerState::OK;
		case 'E':
			return WorkerState::NO_WORK;
		default:
			return WorkerState::FAIL;
		}
	}

	// terminate unresponsive worker
	auto now = std::chrono::steady_clock::now();
	if ((now - m_last_heard) > WORKER_TIMEOUT) {
		SPLOG("Killing unresponsive worker process, pid: %d", m_pid);
		terminate();
		return WorkerState::FAIL;
	}

	if (last_char == 'U') {
		return WorkerState::UPDATE;
	}
	return RUNNING;
}

int WorkerProcess::wait(bool hang)
{
	int options = 0;
	if (!hang) {
		options = WNOHANG;
	}
	int status = 0;
	int retcode;
	while (true) {
		retcode = ::waitpid(m_pid, &status, options);
		if (retcode != -1 || errno != EINTR) {
			break;
		}
	}
	return retcode;
}

void WorkerProcess::terminate()
{
	::kill(m_pid, SIGKILL);
	wait(true);
}

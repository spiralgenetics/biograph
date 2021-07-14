#include "modules/main/main.h"
#include "modules/pipeline/primitives.h"
#include "modules/web/httpclient.h"
#include "modules/mapred/task_attempt.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/task_worker.h"
#include "modules/mapred/task_runner.h"
#include "modules/io/log.h"
#include "modules/io/file_io.h"
#include "modules/io/config.h"

extern "C"
{
	typedef void (*vvfunc)();
};

update_task_runner* g_runner;

void do_notify()
{
	try
	{
		g_runner->update_progress(0.0);
	}
	catch (...)
	{
	}
}

class WorkerMain : public Main
{
public:
	WorkerMain()
	{
		m_usage =
			"Usage:\n"
			"    %s [options] [profile]\n"
		;
	}

protected:
	void add_args() override
	{
		m_options.add_options()
			("profile", po::value(&m_profile), "specify a profile to use")
		;
		m_positional.add("profile", 1);
	}

	int run(po::variables_map vars) override
	{
		return do_worker(m_profile);
	}

private:
	std::string m_profile;
};

std::unique_ptr<Main> worker_main()
{
	return std::unique_ptr<Main>(new WorkerMain);
}

int do_worker(const std::string& profile)
{
	add_primitives();

	int retries = 5;
	// Himem tasks only get one try
	if(profile == "himem") {
		retries = 0;
	}
	task_worker tw(new_taskdb_couch(), retries);

	task_attempt ta;
	if (!tw.get_attempt_for_profile(ta, profile)) {
		printf("E");  // Empty
		return 0;
	}

	SPLOG_P(LOG_DEBUG, "anchored/do_worker> Running '%s' task '%s' with state path %s", profile.c_str(), ta.type.c_str(), ta.state_path.url().c_str());

	task_attempt_result tar;
	update_task_runner runner(tw, ta, tar);
	g_runner = &runner;
	runner.run();

	SPLOG_P(LOG_DEBUG, "anchored/do_worker> Write result = %d", tar.result);
	size_t delay = 2;
	for (size_t i = 0; i < 5; i++) {
		try {
			tw.apply_results(tar);
			break;
		}
		catch(io_exception& io) {
		}
		sleep(delay);
		delay *= 2;
	}

	printf("S");  // Self

	// SPLOG("anchored/do_worker> done");
	return 0;
}

#include "modules/main/main.h"
#include "modules/pipeline/primitives.h"
#include "modules/mapred/taskdb.h"
#include "modules/io/config.h"
#include "modules/io/log.h"

#include <chrono>
#include <csignal>

static
void handle_signals()
{
	for (auto signal : { SIGTERM, SIGHUP }) {
		std::signal(signal, [](int signal) {
			SPLOG("Caught signal: %d", signal);
			http_server::get().stop();
			taskdb_stop_persister();
			std::exit(0);
		});
	}
}

class TaskdbMain : public Main
{
protected:
	int run(po::variables_map vars) override;
};

int TaskdbMain::run(po::variables_map vars)
{
	try {
		SPLOG("Starting taskdb");

		add_primitives();

		taskdb taskdb;
		try {
			taskdb.restore_global_state();
		}
		catch (const io_exception& io) {
			SPLOG("Error when trying to restore taskdb: %s", io.message().c_str());
		}

		// this will periodically save the global state or save when a termination signal is caught.
		int period = CONF_T(int, taskdb_backup_period_in_seconds);
		if (period) {
			taskdb_start_persister(taskdb, std::chrono::seconds(period));
		}

		handle_signals();

		taskdb.register_handlers();

		run_restful_server(CONF(taskdb_bind_list), "", "", "thread");

		return 0;
	}
	catch (const io_exception& io) {
		http_server::get().stop();
		throw;
	}
}

std::unique_ptr<Main> taskdb_main()
{
	return std::unique_ptr<Main>(new TaskdbMain);
}

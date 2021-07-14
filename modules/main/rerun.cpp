#include "modules/main/main.h"
#include "modules/pipeline/primitives.h"
#include "modules/mapred/task_attempt.h"
#include "modules/io/config.h"
#include "modules/io/log.h"

class RerunMain : public Main
{
public:
	RerunMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% --type [task] --path [task path]\n\n"
			"Rerun specific task.\n"
		;
	}
protected:
	void add_args() override;
	int run(po::variables_map vars) override;
	bool needs_cleanup() override { return false; }

private:
	std::string m_task_type;
	std::string m_task_path;
};

void RerunMain::add_args()
{
	m_options.add_options()
		("type", po::value(&m_task_type), "specify the type of the task")
		("path", po::value(&m_task_path), "specify the path of the task")
	;
}

int RerunMain::run(po::variables_map vars)
{
	initialize_app("/reference/");
	log_init("rerun", 2);
	log_build_stamp();
	add_primitives();

	task_attempt ta;
	ta.working_path = CONF_S(path_bulkdata);
	ta.task_id = "test";
	ta.state_counter = 0;
	ta.attempt = 0;
	ta.working_path = path("/tmp/").append_unique("rerun");
	ta.type = m_task_type;
	ta.state_path = m_task_path;

	SPLOG("RerunMain::run> Running task of type: %s", ta.type.c_str());
	try {
		task_attempt_result tar = attempt_task(ta);
		SPLOG("RerunMain::run> Write result = %d", tar.result);
		SPLOG("RerunMain::run> Details: %s", json_serialize(tar).c_str());
	} catch(const io_exception& e) {
		SPLOG("RerunMain::run> Caught io_exception: %s", e.what());
	} catch(...) {
		SPLOG("RerunMain::run> unexpected exception type");
	}

	SPLOG("RerunMain::run> Done");

	return 0;
}

std::unique_ptr<Main> rerun_main()
{
	return std::unique_ptr<Main>(new RerunMain);
}

#include "modules/main/main.h"
#include "modules/pipeline/primitives.h"
#include "modules/mapred/task_attempt.h"
#include "modules/io/config.h"
#include "modules/io/log.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/json_transfer.h"

class DumpTDBMain : public Main
{
public:
	DumpTDBMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% --tmp /path/to/tempdir/\n\n"
			"Dump all of taskdb to stdout.\n"
		;
	}
protected:
	std::string m_parent;
	bool m_full;
	void add_args() override;
	int run(po::variables_map vars) override;
	bool needs_cleanup() override { return false; }
};

void DumpTDBMain::add_args()
{
	m_options.add_options()
		("parent", po::value(&m_parent)->default_value(""), "Parent task to dump subtasks of")
		("full", po::bool_switch(&m_full)->default_value(false), "Dump full json records")
	;

}

int DumpTDBMain::run(po::variables_map vars)
{
	if (m_tmp_dir == "" || path(m_tmp_dir).exists() != path::e_directory) {
		throw io_exception("dump_taskdb needs an existing tmp directory to use");
	}

	initialize_app("");
	auto filename = taskdb_backup_filename(CONF_S(storage_root));
	SPLOG("DumpTDBMain::run> Restoring global state from %s", filename.c_str());
	path backup(filename);
	if (backup.exists() != path::e_file) {
		throw io_exception("taskdb::restore_global_state> Taskdb backup not found");
	}

	task_map_t tasks;
	msgpack_deserialize(tasks, backup.get());
	printf("[\n");
	bool prev = false;
	for(const auto& kvp : tasks) {
		const task_info& task = kvp.second;
		if (task._id.find(m_parent) != 0) {
			continue;
		}
		if (!m_full && task.parent_id != m_parent) {
			continue;
		}
		if (prev) { printf(",\n"); }
		printf("    ");
		if (m_full) {
			printf("%s", json_serialize(task).c_str());
		} else {
			printf("\"_id\":\"%s\", \"type\":\"%s\", \"subtype\":\"%s\", \"state_path\":\"%s\", \"output_path\":\"%s\"",
				task._id.c_str(), task.type.c_str(), task.subtype.c_str(),
				task.state_path.bare_path().c_str(), task.output_path.bare_path().c_str());
		}
		prev = true;
	}
	printf("\n]\n");
	return 0;
}

std::unique_ptr<Main> dump_taskdb_main()
{
	return std::unique_ptr<Main>(new DumpTDBMain);
}

#include "modules/main/main.h"
#include "modules/pipeline/primitives.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/task_attempt.h"
#include "modules/bio_format/exporter.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/json_transfer.h"
#include "modules/bio_format/assembly.h"
#include "modules/bio_format/vcf.h"
#include <fstream>

class ExportMain : public Main
{
public:
	ExportMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1%\n"
			"Given a task from task_db, exports the n'th manifest as type into output.\n"
		;
	}
protected:
	std::string m_task;
	std::string m_type;
	std::string m_out;
	std::string m_param;
	int m_which;
	void add_args() override;
	int run(po::variables_map vars) override;
	bool needs_cleanup() override { return false; }
};

void ExportMain::add_args()
{
	m_options.add_options()
		("task", po::value(&m_task)->required(), "Task to dump output of")
		("type", po::value(&m_type)->required(), "Type of exporter to use")
		("out", po::value(&m_out)->required(), "Output filename")
		("which", po::value(&m_which)->default_value(-1), "Which manifest to dump (-1 means not a vector)")
		("param", po::value(&m_param)->default_value(""), "Parameter to pass exporter")
	;

}

int ExportMain::run(po::variables_map vars)
{
	if (m_tmp_dir == "" || path(m_tmp_dir).exists() != path::e_directory) {
		throw io_exception("export needs an existing tmp directory to use");
	}
	// Get reference directory
        json_spirit::mValue value;
        std::string cfg_file_name = m_tmp_dir + "/config.json";
        std::ifstream cfg_file(cfg_file_name);
        try {
                json_spirit::read_or_throw(cfg_file, value);
        }
        catch(json_spirit::Error_position e) {
		throw io_exception("Couldn't read config file from tmp directory");
        }
        cfg_file.close();
	std::string ref_dir = "";
	try {
       		json_spirit::mObject obj = value.get_obj();
		ref_dir = obj["refdir"].get_str();
	} catch(const std::exception& e) {
		throw io_exception("Unable to recover reference directory");
	}

	fprintf(stderr, "Using reference directory %s\n", ref_dir.c_str());
	initialize_app(ref_dir);
	auto filename = taskdb_backup_filename(CONF_S(storage_root));
	path backup(filename);
	if (backup.exists() != path::e_file) {
		throw io_exception("ExportMain::run> Taskdb backup not found");
	}

	task_map_t tasks;
	msgpack_deserialize(tasks, backup.get());
	const auto& it = tasks.find(m_task);
	if (it == tasks.end()) {
		throw io_exception(printstring("ExportMain::run> No such task %s", m_task.c_str()));
	}
	std::string data = it->second.output_path.get();
	manifest m;
	try {
		json_deserialize(m, data);
	} catch (const io_exception& e) {
		std::vector<manifest> mv;
		try {
			json_deserialize(mv, data);
		}
		catch(const io_exception& e) {
			throw io_exception("Task output is not a manifest or std::vector<manifest>");
		}
		if (m_which == -1) {
			m_which = 0;
			if (mv.size() != 1) {
				throw io_exception("Task output has multiple manifests and which not specified");
			}
		}
		if (m_which < 0 || m_which >= int(mv.size())) {
			throw io_exception("Which is out of range");
		}
		m = mv[m_which];
	}
	fprintf(stderr, "Reading from manifest, # of records = %lu, # of bytes = %lu\n",
		m.get_num_records(), m.get_size());

	manifest_reader mr(m);
	std::unique_ptr<exporter> e;
	file_writer fw(m_out);

        if (m_type == "vcf_bp") {
		e = make_unique<vcf_exporter>(fw, "", false);
        } else {
		fprintf(stderr, "Making exporter\n");
		e = registry<exporter, writable&, bool, const std::string&>::get_safe(m_type, fw, false, m_param);
		fprintf(stderr, "Calling exporter\n");
        }
	e->export_from(mr);

	return 0;
}

std::unique_ptr<Main> export_main()
{
	return std::unique_ptr<Main>(new ExportMain);
}

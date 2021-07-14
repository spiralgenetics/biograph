#pragma once

#include <map>
#include <memory>
#include <string>
#include <boost/program_options.hpp>
#include "modules/io/version.h"
#include "modules/io/runtime_stats.h"
#include "modules/mapred/taskdb.h"

namespace po = boost::program_options;

class Main;
typedef std::function<std::unique_ptr<Main>()> main_f;

class Main
{
public:
	int main(const std::string& name, int argc, char* argv[]);
  virtual ~Main() = default;

	// default version for worker, resurrect, etc. misc commands
	virtual const product_version& get_version()  { return biograph_current_version; }
	virtual po::variables_map parse_args(int argc, char* argv[]);
	virtual int run(po::variables_map vars) { return 0; }

	virtual void print_help(std::ostream& os, bool show_advanced = false);
	virtual void add_args() {}
	void cleanup(bool success = true);
	virtual bool needs_cleanup() { return true; }

protected:
	pid_t launch(const char* name, main_f fn, std::vector<const char*> args = {});
    int main_child(const char* name, main_f fn, std::vector<const char*> args = {});
	void initialize_app(const std::string& ref_dir, const std::string& log_file = "");
	void launch_daemons();

	std::string m_name;
	std::string m_usage =
		"Usage:\n"
		"    %s [options]\n"
	;

	// Number of requested threads, or "auto"
	std::string m_requested_threads;
	// Number of actual threads
	size_t m_num_threads;

	std::string m_tmp_dir;
	bool m_tmp_dir_made = false;
	bool m_keep_tmp;
	bool m_debug_log = false;
	bool m_cache_all = false;
	std::string m_log_file;
	std::string m_stats_file;
#if GPERFTOOLS
	std::string m_cpuprofile_dir;
#endif
	runtime_stats m_stats;

	unsigned int columns = get_terminal_width();

	po::positional_options_description m_positional;

	po::options_description m_general_options{"General Options", columns};
	po::options_description m_kmer_options{"Kmerization Options", columns};
	po::options_description m_correction_options{"Read Correction Options", columns};
	po::options_description m_variant_options{"Variant Calling Options", columns};
	po::options_description m_assembly_options{"Assembly Options", columns};
	po::options_description m_advanced_options{"Advanced Options", columns};
	po::options_description m_secret_options{"Not included in help", columns};

	po::options_description m_all_options{"All Options", columns};
	po::options_description m_options{"", columns};

	void set_mem_limit(uint64_t max_mem);
	size_t get_free_port();

	std::string m_cmdline;

private:
	pid_t m_logger_pid = 0;
	pid_t m_taskdb_pid = 0;
	pid_t m_normal_pid = 0;
	pid_t m_himem_pid = 0;
};

std::unique_ptr<Main> rerun_main();
std::unique_ptr<Main> dump_taskdb_main();
std::unique_ptr<Main> resurrect_main();
std::unique_ptr<Main> export_main();

int do_worker(const std::string& profile);
int do_migration(const std::string& root_url, bool fake_mode);

#include <fstream>
#include <signal.h>
#include <stdexcept>
#include <sys/prctl.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>

#include "modules/bio_base/dna_base_set.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_format/fastq.h"
#include "modules/bio_format/vcf.h"
#include "modules/bio_mapred/compute_coverage.h"
#include "modules/bio_mapred/kmerize_bf.h"
#include "modules/bio_mapred/mem_seqset.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/utils.h"
#include "modules/io/zip.h"
#include "modules/mapred/dual_map_task.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/task_mgr.h"
#include "modules/pipeline/dataset_path.h"
#include "modules/pipeline/paired_merger.h"

#include "modules/main/main.h"

using namespace boost::algorithm;

class RefToSEQSETMain: public Main
{
	using contig_remap_t = std::map<int, int>;

public:
	RefToSEQSETMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% [OPTIONS] --ref [ref dir] --out [seqset]\n\n"
			"Generate a BioGraph seqset based on a reference."
		;
	}

protected:
	int run(po::variables_map vars) override;
	void add_args() override;

	template<class OutType>
	void run_task(OutType& out, std::unique_ptr<task> task);

	void check_for_terminate();

	const product_version& get_version() override { return biograph_current_version; }

private:
	std::string m_read_size;
	std::string m_out_file;
	std::string m_ref_dir;
	bool m_force;
};

static volatile bool terminate = false;

// call this to see if it's time to terminate
void RefToSEQSETMain::check_for_terminate()
{
	if(terminate) {
		std::cout << "\nControl-C detected.\n";
		SPLOG("Control-C detected.");
		m_keep_tmp = true;
		cleanup(false);
		::exit(1);
	}
}

// anchored handles termination in the main loop.
static void signal_handler(int sig)
{
	// One is enough
	signal(sig, SIG_IGN);
	terminate = true;
}

// helper to only update progress when the delta is > 0.01%
static void update_progress(const float& cur_progress, float& prev_progress) {
	if(cur_progress - prev_progress > 0.0001) {
		prev_progress = cur_progress;
		print_progress(cur_progress);
	}
}

void RefToSEQSETMain::add_args()
{
	m_general_options.add_options()
		("ref", po::value(&m_ref_dir)->default_value(""), "Reference directory")
		("read-size", po::value(&m_read_size)->required(), "Read size for seqset")
		("out", po::value(&m_out_file)->required(), "Output seqset file")
		("force,f", po::bool_switch(&m_force)->default_value(false), "Overwrite existing output file")
	;
    m_general_options.add(track_mem_program_options());

	m_options.add(m_general_options);
}

static size_t validate_param(const std::string& param, const std::string& value, const std::vector<size_t>range = {})
{
	size_t num_val;
	try {
		num_val = std::stoull(value);
	}
	catch (const std::exception& ex) {
		throw std::runtime_error(boost::str(boost::format("%s must specify an integer") % param));
	}
	if((range.size() >= 1) && (num_val < range[0])) {
		throw std::runtime_error(boost::str(boost::format("%s must specify an integer >= %d") % param % range[0]));
	}
	if((range.size() == 2) && (num_val > range[1])) {
		throw std::runtime_error(boost::str(boost::format("%s must specify an integer <= %d") % param % range[1]));
	}
	return num_val;
}

template<class OutType>
void RefToSEQSETMain::run_task(OutType& out, std::unique_ptr<task> task)
{
	task_mgr tm(new_taskdb_couch());
	std::string id = tm.add_job(CONF_S(path_bulkdata), std::move(task), "ref2seqset");
	int job_state = 0;
	int tdb_errs = 0;
	float prev_progress = 0.0;
	float cur_progress = 0.0;

	print_progress(0.0);
	while(job_state == 0) {
		try {
			job_state = tm.state(id);
			cur_progress = tm.get_progress(id);
			update_progress(cur_progress, prev_progress);
			tdb_errs = 0;
		} catch(...) {
			tdb_errs++;
		}
		if (tdb_errs > 5) {
			throw io_exception("Can't communicate with taskdb");
		}

		usleep(500000);
		check_for_terminate();
	}
	if (job_state != 1) {
		throw io_exception(tm.get_error(id));
	}
	print_progress(1.0);
	std::cout << "\n";
	tm.get_output(out, id);
}

int RefToSEQSETMain::run(po::variables_map vars)
{
	if(boost::filesystem::exists(m_out_file) && !m_force) {
		std::cerr << "Refusing to overwrite '" << m_out_file << "'. Use -f to override.\n";
		exit(1);
	}

	size_t read_size = validate_param("--read-size", m_read_size, {30, 255});

	// Initialize and kick off the daemons
	initialize_app(m_ref_dir);
	launch_daemons();

	// Now set up the custom handler
	signal(SIGINT, signal_handler);
	signal(SIGTERM, signal_handler);

	std::cout << "Running seqset generation\n";
	auto seqset_task = make_unique<mem_seqset_task>();
	seqset_task->read_size = read_size;
	seqset_task->num_threads = m_num_threads;
	seqset_task->max_mem = get_maximum_mem_bytes() / 1024 / 1024;
	manifest seqset;
	run_task(seqset, std::move(seqset_task));

	path out = *seqset.begin();
	boost::filesystem::copy_file(out.bare_path(), m_out_file,
		boost::filesystem::copy_option::overwrite_if_exists);
	return 0;
}

std::unique_ptr<Main> ref2seqset_main()
{
	return std::unique_ptr<Main>(new RefToSEQSETMain);
}



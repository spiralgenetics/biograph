#include <boost/filesystem.hpp>
#include <stdexcept>

#include "modules/main/main.h"

#include "modules/bio_base/bwt_file.h"
#include "modules/bio_mapred/make_bwt.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_buffer.h"
#include "modules/mapred/task_mgr.h"

class RefToBWTMain: public Main
{
public:
	RefToBWTMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% [OPTIONS] --in [file.ref] --out [file.bwt]\n\n"
			"Convert a .ref to .bwt\n"
		;
	}

protected:
	const product_version& get_version() override { return biograph_current_version; }
	void add_args() override;
	int run(po::variables_map vars) override;

private:
	std::string m_ref_file;
	std::string m_bwt_file;
	size_t m_cent_mod;
};

void RefToBWTMain::add_args()
{
	m_options.add_options()
		("in", po::value(&m_ref_file)->required(), "The input ref file")
		("out", po::value(&m_bwt_file)->required(), "The output bwt file")
		("cent_mod", po::value(&m_cent_mod)->default_value(64), "Make a century table entry every cent_mod records")
	;
	m_positional.add("in", 1);
	m_positional.add("out", 1);
}

int RefToBWTMain::run(po::variables_map vars)
{
	initialize_app("");
	launch_daemons();

	task_mgr tm(new_taskdb_couch());

	auto bwt_task = make_unique<make_bwt_task>();
	bwt_task->input_ref = m_ref_file;
	bwt_task->output_bwt = m_bwt_file;
	bwt_task->cent_mod = m_cent_mod;

	std::string id = tm.add_job(CONF_S(path_bulkdata), std::move(bwt_task), "ref2bwt");

	std::string out;

	int job_state = 0;
	int tdb_errs = 0;
	while(job_state == 0) {
		try {
			job_state = tm.state(id);
			// cur_progress = tm.get_progress(id);
			// update_progress(cur_progress, prev_progress);
			tdb_errs = 0;
		} catch(...) {
			tdb_errs++;
		}
		if (tdb_errs > 5) {
			throw io_exception("Can't communicate with taskdb");
		}

		usleep(500000);
		// check_for_terminate();
	}
	if (job_state != 1) {
		throw io_exception(tm.get_error(id));
	}


	tm.get_output(out, id);

	std::cout << out << " saved.\n";

	return 0;
}

std::unique_ptr<Main> ref2bwt_main()
{
	return std::unique_ptr<Main>(new RefToBWTMain);
}

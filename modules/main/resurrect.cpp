
#include "modules/main/main.h"
#include "modules/mapred/task_mgr.h"
#include "modules/io/config.h"
#include "modules/io/log.h"

class ResurrectMain : public Main
{
public:
	ResurrectMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% --refdir [ref_dir] --job [job_id] --tmp [tmp_dir]\n\n"
			"Resurrect a job.\n"
		;
	}
protected:
	void add_args() override;
	int run(po::variables_map vars) override;

private:
	std::string m_ref_dir;
	std::string m_job_id;
};

void ResurrectMain::add_args()
{
	m_options.add_options()
		("refdir", po::value(&m_ref_dir)->default_value(""), "Reference directory created by make_ref")
		("job", po::value(&m_job_id), "Job id to resurrect")
	;
}

int ResurrectMain::run(po::variables_map vars)
{
        if(m_ref_dir.empty() || m_job_id.empty() || m_tmp_dir.empty()) {
                print_help(std::cerr, true);
                ::exit(1);
        }

	initialize_app(m_ref_dir);
	launch_daemons();
        task_mgr tm(new_taskdb_couch());
        tm.resurrect_job(m_job_id);
        int job_state = 0;
        int tdb_errs = 0;
        while(job_state == 0) {
                try {
                        job_state = tm.state(m_job_id);
                        double progress = tm.get_progress(m_job_id);
                        printf("Progress = %f\n", progress);
                        tdb_errs = 0;
                } catch(...) {
                        tdb_errs++;
                }
                if (tdb_errs > 5) {
                        throw io_exception("Can't communicate the taskdb");
                }
                sleep(1);
        }
        if (job_state != 1) {
                throw io_exception(tm.get_error(m_job_id));
        }
	printf("Job completed succesfully\n");

	return 0;
}

std::unique_ptr<Main> resurrect_main()
{
	return std::unique_ptr<Main>(new ResurrectMain);
}

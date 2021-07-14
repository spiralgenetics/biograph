#include <signal.h>
#include <stdexcept>
#include <sys/prctl.h>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/version.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_mapred/make_readmap.h"

#include "modules/main/main.h"


class MigrateReadmapMain: public Main
{
public:
	MigrateReadmapMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% [OPTIONS] --orig-seqset [seqset] --orig-readmap [readmap] --new-seqset [seqset] --new-readmap [readmap]\n\n"
			"Migrates a readmap from one sequence set (orig) to another (new)\n";
		;
	}

protected:
	void add_args() override;
	int run(po::variables_map vars) override;
	const product_version& get_version() override { return biograph_current_version; }

private:
	std::string m_orig_seqset;
	std::string m_orig_readmap;
	std::string m_new_seqset;
	std::string m_new_readmap;
	bool m_force;
};

// anchored handles termination in the main loop.
static void signal_handler(int sig)
{
	// One is enough
	signal(sig, SIG_IGN);
	std::cout << "\nControl-C detected.\n";
	::exit(1);
}


void MigrateReadmapMain::add_args()
{
	m_general_options.add_options()
		("orig-seqset", po::value(&m_orig_seqset)->required(), "Original sequence set for original readmap")
		("orig-readmap", po::value(&m_orig_readmap)->required(), "Original readmap")
		("new-seqset", po::value(&m_new_seqset)->required(), "New sequence set to migrate to")
		("new-readmap", po::value(&m_new_readmap)->required(), "New output readmap")
		("force,f", po::bool_switch(&m_force)->default_value(false), "Overwrite existing output file")
	;
	m_options.add(m_general_options);
}

int MigrateReadmapMain::run(po::variables_map vars)
{
	if(boost::filesystem::exists(m_new_readmap)) {
		if (m_force) {
			boost::filesystem::remove(m_new_readmap);
		} else {
			std::cerr << "Refusing to overwrite '" << m_new_readmap << "'. Use -f to override.\n";
			exit(1);
		}
	}

	initialize_app("");

	// initialize_app() ignores SIGINT, so handle it ourselves.
	signal(SIGINT, signal_handler);

	std::cerr << "Loading original sequence set\n";
    auto orig_seqset = std::make_shared<seqset>(m_orig_seqset);
	std::cerr << "Loading new sequence set\n";
    auto new_seqset = std::make_shared<seqset>(m_new_seqset);
	std::cerr << "Loading original readmap\n";
    readmap orig_readmap(new_seqset, m_orig_readmap);
	std::cerr << "Doing migration\n";
    make_readmap::migrate(*orig_seqset, orig_readmap, *new_seqset, m_new_readmap, false);

	return 0;
}

std::unique_ptr<Main> migrate_readmap_main()
{
	return std::unique_ptr<Main>(new MigrateReadmapMain);
}



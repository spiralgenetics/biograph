#include <boost/filesystem.hpp>
#include <stdexcept>

#include "modules/main/main.h"

#include "modules/bio_base/bwt_file.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_buffer.h"

class BWTQueryMain: public Main
{
public:
	BWTQueryMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% [OPTIONS] --in [file.bwt] --query [DNA String]\n\n"
			"Query a BWT reference for a given kmer.\n"
		;
	}

protected:
	const product_version& get_version() override { return biograph_current_version; }
	void add_args() override;
	int run(po::variables_map vars) override;

private:
	std::string m_bwt_file;
	std::string m_query_kmer;
	bool m_verbose;
	bool m_quiet;

	void query_bwt();
};

void BWTQueryMain::add_args()
{
	m_options.add_options()
		("in", po::value(&m_bwt_file)->required(), "Reference bwt file to search")
		("query", po::value(&m_query_kmer)->required(), "Query kmer, e.g. \"AGTTCGA\"")
		("verbose", po::bool_switch(&m_verbose)->default_value(false), "Output more than 10 prefixes"
			" (could produce large outputs!)")
		("quiet", po::bool_switch(&m_quiet)->default_value(false), "Only output the graph kmers and warnings or errors")
	;
}

int BWTQueryMain::run(po::variables_map vars)
{
	query_bwt();
	return 0;
}

void BWTQueryMain::query_bwt()
{
	bwt_file file(m_bwt_file);
	bwt_range query_range = file.bwt().find(dna_sequence(m_query_kmer));

	if (! m_quiet) std::cerr << "Query: \"" << m_query_kmer << "\"\n";

	if (! query_range.valid()) {
		if (! m_quiet) std::cerr << "No valid results were found.\n";
		return;
	}

	uint64_t entry_count = query_range.matches();
	if (! m_quiet) std::cerr << "Found " << entry_count << " entries\n";

	if (entry_count > 10 && (! m_verbose)) {
		std::cerr << "More than ten entries matched the query.  Use the \"--verbose\" option to see them all.\n";
		entry_count = 10;
	}

	for(size_t i = 0; i < entry_count; i++) {
		std::cout << query_range.get_match(i) << "\n";
	}
}

std::unique_ptr<Main> bwt_query_main()
{
	return std::unique_ptr<Main>(new BWTQueryMain);
}

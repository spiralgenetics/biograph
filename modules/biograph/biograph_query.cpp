
#include "modules/main/main.h"

#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_buffer.h"
#include "modules/bio_format/vcf.h"
#include "modules/bio_base/seqset.h"
#include <stdexcept>

#include <boost/filesystem.hpp>

class SEQSETQueryMain: public Main
{
public:
	SEQSETQueryMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% [OPTIONS] --in [file.seqset] --query [DNA String]\n\n"
			"Query a seqset for a given kmer.\n"
		;
	}

protected:
	const product_version& get_version() override { return biograph_current_version; }
	void add_args() override;
	int run(po::variables_map vars) override;

private:
	std::string m_seqset_file;
	std::string m_query_kmer;
	bool m_verbose;
	bool m_quiet;

	void query_seqset();
};

void SEQSETQueryMain::add_args()
{
	m_options.add_options()
		("in", po::value(&m_seqset_file)->required(), "Get the prefixes from this seqset file")
		("query", po::value(&m_query_kmer)->required(), "Query kmer, e.g. \"AGTTCGA\"")
		("verbose", po::bool_switch(&m_verbose)->default_value(false), "Output more than 10 prefixes"
			" (could produce large outputs!)")
		("quiet", po::bool_switch(&m_quiet)->default_value(false), "Only output the graph kmers and warnings or errors")
	;
}

int SEQSETQueryMain::run(po::variables_map vars)
{
	query_seqset();
	return 0;
}

void SEQSETQueryMain::query_seqset()
{
	seqset_file file(m_seqset_file);
	const seqset& the_seqset = file.get_seqset();
	seqset_range query_context = the_seqset.find(m_query_kmer);
	if (! m_quiet) std::cerr << "Query: \"" << m_query_kmer << "\"\n";

	if (! query_context.valid()) {
		if (! m_quiet) std::cerr << "No valid results were found.\n";
		return;
	}

	if (! m_quiet) std::cerr << "Query has " << static_cast<int>(query_context.size()) << " bases.\n";

	uint64_t entry_count = query_context.end() - query_context.begin();
	if (! m_quiet) std::cerr << "Found " << entry_count << " entries\n";

	if (entry_count > 10 && (! m_verbose)) {
		std::cerr << "More than ten entries matched the query.  Use the \"--verbose\" option to see them all.\n";
		entry_count = 10;
	}
    uint64_t seqset_id = query_context.begin();
	while(entry_count) {
      query_context = the_seqset.ctx_entry(seqset_id);
      std::cout << query_context.sequence().as_string() << "\n";
      ++seqset_id;
      entry_count--;
	}
}

std::unique_ptr<Main> seqset_query_main()
{
	return std::unique_ptr<Main>(new SEQSETQueryMain);
}

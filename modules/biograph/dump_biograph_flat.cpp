
#include "modules/main/main.h"

#include "modules/io/log.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_format/dna_io.h"

#include <stdexcept>


class SEQSETDumpFlatMain: public Main
{
public:
	SEQSETDumpFlatMain()
	{
		m_usage =
			"%1% version %2%\n\n"
			"Usage: %1% [OPTIONS] --in [file.seqset.flat]\n\n"
			"Dump to stdout a string representation ofa flat file\n"
		;
	}

protected:
	const product_version& get_version() override { return biograph_current_version; }
	void add_args() override;
	int run(po::variables_map vars) override;

private:
	std::string m_flat_file;

	void seqset_dump_flat();
};

void SEQSETDumpFlatMain::add_args()
{
	m_options.add_options()
		("in", po::value(&m_flat_file)->required(), "Flat file to get seqs from")
	;
}

int SEQSETDumpFlatMain::run(po::variables_map vars)
{
	seqset_dump_flat();

	return 0;
}

void SEQSETDumpFlatMain::seqset_dump_flat()
{

	dna_reader din (make_unique<file_reader>(m_flat_file));
	dna_sequence seq = din.read();
	while (seq.size()) {
		std::cout << seq.as_string() << std::endl; 
		seq = din.read();
	}
}

std::unique_ptr<Main> seqset_dump_main()
{
	return std::unique_ptr<Main>(new SEQSETDumpFlatMain);
}

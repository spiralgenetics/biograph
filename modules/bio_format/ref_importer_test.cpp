#include "modules/test/test_utils.h"
#include "modules/bio_format/fasta_ref_importer.h"
#include "gtest/gtest.h"

TEST(reference_importer, e_coli)
{
	std::vector<std::string> nothing;
	file_reader e_coli("golden/e_coli.fasta");
	fasta_ref_importer ri(path(make_path("e_coli_ref/")), e_coli, nothing, 50, no_update);
	ri.run();
}

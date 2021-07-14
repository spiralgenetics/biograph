#include "gtest/gtest.h"
#include "modules/pipeline/primitives.h"
#include "modules/bio_format/fasta_ref_importer.h"
#include "modules/bio_base/reference.h"
#include "modules/test/test_utils.h"
#include "modules/test/build_ref.h"

TEST(reference_wrapper_test, basic)
{
	using flat_ref_range_t = std::pair<size_t, size_t>;

	std::vector<std::string> nothing;
	file_reader e_coli("golden/fake_ref.fasta");
	std::string fake_ref_dir_path = make_path("fake_ref/");
	fasta_ref_importer ri(path(fake_ref_dir_path), e_coli, nothing, 50, no_update);
	ri.run();
	perform_build_ref("fake", "golden/fake_ref.fasta");

	reference fake_ref("fake");

	flat_ref_range_t the_As = fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 60, 120, true);
	dna_const_iterator seq_start = fake_ref.get_dna(the_As.first);
	dna_const_iterator seq_end = fake_ref.get_dna(the_As.second);

	std::string sixty_As = dna_sequence(seq_start, seq_end).as_string();
	ASSERT_EQ(sixty_As, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

	ASSERT_THROW(fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 59, 120, true), io_exception);
	ASSERT_THROW(fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 60, 121, true), io_exception);
	ASSERT_THROW(fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 59, 121, true), io_exception);

	ASSERT_THROW(fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 119, 110, true), io_exception);

	the_As = fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 59, 120, false);
	seq_start = fake_ref.get_dna(the_As.first);
	seq_end = fake_ref.get_dna(the_As.second);
	sixty_As = dna_sequence(seq_start, seq_end).as_string();
	ASSERT_EQ(sixty_As, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

	the_As = fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 60, 121, false);
	seq_start = fake_ref.get_dna(the_As.first);
	seq_end = fake_ref.get_dna(the_As.second);
	sixty_As = dna_sequence(seq_start, seq_end).as_string();
	ASSERT_EQ(sixty_As, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

	ASSERT_THROW(fake_ref.flatten_range("Sixty_Bases_Separated_by_60_Ns", 59, 121, false), io_exception);
}

#include <gtest/gtest.h>

#include "modules/test/test_utils.h"

#include "modules/bio_base/flat_ref.h"
#include "modules/mapred/path.h"


// Class to expose and gut the flat_ref for testing purposes.
class flat_ref_tester
{
	flat_ref& m_flat_ref_ref;

public:
	flat_ref_tester(flat_ref& flat_ref_ref)
		: m_flat_ref_ref(flat_ref_ref) {}

	// Caution!  This function guts the original flat_ref!
	std::unique_ptr<flat_ref::index_t> expose_index() const
	{
		return std::move(m_flat_ref_ref.m_index);
	}

	mmap_buffer&& get_dna_mmap() { return std::move(m_flat_ref_ref.m_mmap); }
};


TEST(flat_ref, bad)
{
	std::string golden_fasta_path("golden/bad.fasta");
	std::string flat_ref_path(make_path("bad.flat"));

	file_reader fasta_reader(golden_fasta_path);
	file_writer flat_writer(flat_ref_path);
	SPLOG("About to build flat ref.");
	flat_ref_builder frb(fasta_reader, flat_writer);
	SPLOG("Flat ref builder constructed.");
	ASSERT_THROW(frb.run(), io_exception);
}


TEST(flat_ref, round_trip)
{
	std::string golden_fasta_path("golden/sequences.fasta");
	std::string flat_ref_path(make_path("sequences.flat"));

	file_reader fasta_reader(golden_fasta_path);
	file_writer flat_writer(flat_ref_path);
	flat_ref_builder frb(fasta_reader, flat_writer);
	frb.run();

	std::string flattened_fasta_path(make_path("flattened.fasta"));
	flat_ref flattened_sequences(flat_ref_path);
	file_writer flattened_fasta(flattened_fasta_path);
	flattened_sequences.make_fasta(flattened_fasta);
	flattened_fasta.close();

	ASSERT_TRUE(diff(flattened_fasta_path, golden_fasta_path));
}


TEST(flat_ref, IUPAC_round_trip)
{
	std::string golden_fasta_path("golden/random_iupac.fasta");
	std::string flat_ref_path(make_path("sequences.flat"));

	file_reader fasta_reader(golden_fasta_path);
	file_writer flat_writer(flat_ref_path);
	flat_ref_builder frb(fasta_reader, flat_writer);
	frb.run();

	std::string flattened_fasta_path(make_path("random_iupac.fasta"));
	flat_ref flattened_sequences(flat_ref_path);
	file_writer flattened_fasta(flattened_fasta_path);
	flattened_sequences.make_fasta(flattened_fasta, 60);
	flattened_fasta.close();

	ASSERT_TRUE(diff(flattened_fasta_path, golden_fasta_path));
}


TEST(flat_ref, block)
{
	std::string golden_fasta_path("golden/sequences.fasta");
	std::string flat_ref_path(make_path("sequences.flat"));

	file_reader fasta_reader(golden_fasta_path);
	file_writer flat_writer(flat_ref_path);
	flat_ref_builder frb(fasta_reader, flat_writer);
	frb.run();

	std::string flattened_fasta_path(make_path("flattened.fasta"));
	flat_ref flattened_sequences(flat_ref_path);

	flat_ref_tester tester(flattened_sequences);
	mmap_buffer tester_buffer{tester.get_dna_mmap()};
	ASSERT_EQ(std::memcmp(tester_buffer.buffer(), flat_ref::k_magic_header.data(), flat_ref::k_magic_header.size()), 0);
	mem_io block_buffer("", track_alloc("flat_ref_test"));
	block_buffer.write(tester_buffer.buffer() + flat_ref::k_magic_header.size(), tester_buffer.size());
	flat_ref block_constructed(tester.expose_index(), std::move(block_buffer));

	file_writer flattened_fasta(flattened_fasta_path);
	block_constructed.make_fasta(flattened_fasta);
	flattened_fasta.close();

	ASSERT_TRUE(diff(flattened_fasta_path, golden_fasta_path));
}

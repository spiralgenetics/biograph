#include "modules/test/build_ref.h"
#include "modules/pipeline/primitives.h"
#include "modules/bio_mapred/overlap_graph.h"
#include "modules/io/progress.h"
#include "modules/bio_base/reference.h"

#include <vector>
#include <gtest/gtest.h>

struct dup_record 
{
	dup_record(size_t _i, size_t _j, bool _flipped) 
		: i(_i), j(_j), flipped(_flipped) {}
	size_t i;
	size_t j;
	bool flipped;
};

struct overlap_record 
{
	overlap_record(size_t _i, bool _flipped, size_t _overlap) 
		: i(_i), flipped(_flipped), overlap(_overlap) {}
	size_t i;
	bool flipped;
	size_t overlap;
};

TEST(overlap_graph_test, basic) 
{
	std::vector<dna_sequence> seqs;
	seqs.push_back(dna_sequence("ACTACCGTAC")); // 0
	seqs.push_back(dna_sequence("TACCGTACCC")); // 1
	seqs.push_back(dna_sequence("TACCGTACCG")); // 2
	seqs.push_back(dna_sequence("ACCGTACCGT")); // 3
	seqs.push_back(dna_sequence("GTACGGTAGT")); // 4
	seqs.push_back(dna_sequence("CCAATATTGG")); // 5
	seqs.push_back(dna_sequence("AGTACGGTAG")); // 6
	seqs.push_back(dna_sequence("TACTACCGTA")); // 7
	seqs.push_back(dna_sequence("AGTACGGTAG")); // 8

	overlap_graph<std::vector<dna_sequence> > og(seqs);
	noisy_progress_handler nph;
	std::vector<dup_record> dups;
	og.prepare(nph, [&](size_t i, size_t j, bool flipped) {
		SPLOG("Found dup: %zd, %zd, %d", i, j, flipped);
		dups.push_back(dup_record(i, j, flipped));
	});
	ASSERT_EQ(dups.size(), size_t(2));
	ASSERT_EQ(dups[0].i, size_t(0));
	ASSERT_EQ(dups[0].j, size_t(4));
	ASSERT_EQ(dups[0].flipped, true);
	ASSERT_EQ(dups[1].i, size_t(6));
	ASSERT_EQ(dups[1].j, size_t(8));
	ASSERT_EQ(dups[1].flipped, false);

	// Do formal check
	std::vector<overlap_record> overlaps;
	og.find_overlaps(0, true, 6, [&](size_t j, bool flipped, size_t overlap) {
		overlaps.push_back(overlap_record(j, flipped, overlap));
	});
	ASSERT_EQ(overlaps.size(), size_t(4));
	ASSERT_EQ(overlaps[0].i, size_t(6));
	ASSERT_EQ(overlaps[0].flipped, true);
	ASSERT_EQ(overlaps[0].overlap, size_t(9));
	ASSERT_EQ(overlaps[1].i, size_t(1));
	ASSERT_EQ(overlaps[1].flipped, false);
	ASSERT_EQ(overlaps[1].overlap, size_t(8));
	ASSERT_EQ(overlaps[2].i, size_t(2));
	ASSERT_EQ(overlaps[2].flipped, false);
	ASSERT_EQ(overlaps[2].overlap, size_t(8));
	ASSERT_EQ(overlaps[3].i, size_t(3));
	ASSERT_EQ(overlaps[3].flipped, false);
	ASSERT_EQ(overlaps[3].overlap, size_t(7));
	overlaps.clear();
	og.find_overlaps(0, false, 6, [&](size_t j, bool flipped, size_t overlap) {
		overlaps.push_back(overlap_record(j, flipped, overlap));
	});
	ASSERT_EQ(overlaps.size(), size_t(1));
	ASSERT_EQ(overlaps[0].i, size_t(7));
	ASSERT_EQ(overlaps[0].flipped, false);
	ASSERT_EQ(overlaps[0].overlap, size_t(9));

	// Print some values look at for fun.
	for(size_t i = 0; i < 9; i++)
	{
		SPLOG("Checking for forward overlaps for %s", seqs[i].as_string().c_str());
		og.find_overlaps(i, true, 6, [&](size_t j, bool flipped, size_t overlap) {
			SPLOG("  Found overlap: %d:%d %s (%d)", (int) j, flipped, 
				flipped ? seqs[j].rev_comp().as_string().c_str()
				        : seqs[j].as_string().c_str(), 
				(int) overlap);
		});
		SPLOG("Checking for reverse overlaps for %s", seqs[i].as_string().c_str());
		og.find_overlaps(i, false, 6, [&](size_t j, bool flipped, size_t overlap) {
			SPLOG("  Found overlap: %d:%d %s (%d)", (int) j, flipped, 
				flipped ? seqs[j].rev_comp().as_string().c_str()
				        : seqs[j].as_string().c_str(), 
				(int) overlap);
		});
		
	}
}

class ref_as_reads
{
public:
	ref_as_reads(const std::string& refname, size_t read_len) 
		: m_reference(refname)
		, m_read_len(read_len)
	{
		SPLOG("Loading reference");
	}
	size_t size() const { return m_reference.size() - m_read_len; }
	dna_slice operator[](size_t i) const 
	{
		return dna_slice(m_reference.get_dna(i), m_read_len);
	}
private:
	reference m_reference;
	size_t m_read_len;
};

TEST(overlap_graph_test, ref) 
{
	perform_build_ref("e_coli", "datasets/fasta/e_coli_k12.ASM584v1.fasta");
	ref_as_reads rr("e_coli", 100);
	overlap_graph<ref_as_reads> og(rr);
	noisy_progress_handler nph;
	size_t dupcount = 0;
	og.prepare(nph, [&](size_t i, size_t j, bool flipped) {
		dupcount++;
		//SPLOG("Found dup: %zd, %zd, %d", i, j, flipped);
        });
	SPLOG("Found %zd dups", dupcount);
}


#include <gtest/gtest.h>

#include "modules/bio_base/seq_position.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_format/struct_var.h"
#include "modules/bio_base/reference.h"
#include "modules/pipeline/primitives.h"

#include "modules/test/build_ref.h"

TEST(edit_distance, pure_insert)
{
	struct_var	pure_insert;
	pure_insert.is_structural = true;
	pure_insert.ref_start = seq_position(0, 1793333);
	pure_insert.rev_start = false;
	pure_insert.ref_end = seq_position(0, 1793334);
	pure_insert.rev_end = false;
	pure_insert.assembled = dna_sequence {
		"TGCTCTGAAAAGAGTAACGCGCTTTACTATCCCTGACAATCACCAACAACATCGAACAAGAT"
		"AATAAATTCCTGGTTTAATATCCGACAAGTGAAAACATGCACCCGGACGGGCAGCATGTCGCTCCACAAGTGCAGAGCTT"
		"ACTTGTGTTGTACCGAAGCACTCTGTTCAGGTGGCTGATAGTTGTCAATGTGACTCGCCACGCCAAGAAGAATGACTGAA"
		"ACGACAAGAACGATCCAACCTGTTAATTCAATAAGACGATTCATTACAGCCCACATCTCTCTTGATTGATCCATTAACTT"
		"CAGGGGGTAAATGTTACTTAGCAATAATAGCTCCAGCAAGATTTTTACTGAGGTTTTTGCGATATTTAGCTTTTGTCGTT"
		"GGAAAATTCGCTATTTTTTTGACTTAAGTTAAACAACATCCCTTATTGCTGGCAGGTTATTAAACTGTTGAGCGTGGGTA"
		"AAGGATAGTGTCAAATAGCCATCATACTTCAATGAGAGGCAATGACATGAGCGACAACATCCGTGTTGGGTTGATTGGGTATGG"
	};
	pure_insert.var_start = 80;
	pure_insert.var_end = 468;
	pure_insert.depth = 100;
	pure_insert.var_id = 0;
	pure_insert.flipped = false;
	pure_insert.is_ambig = false;
	pure_insert.avg_depth = 100.0;
	pure_insert.min_overlap = 70;
	pure_insert.avg_overlap = 75.2;
	pure_insert.is_ambig = false;
	pure_insert.has_holes = false;
	pure_insert.align_failed = false;
	pure_insert.sub_id = -1;
	
	perform_build_ref("e_coli_dh10b_CP000948.1", "datasets/fasta/e_coli_dh10b_CP000948.1.fasta");
	reference e_coli_ref("e_coli_dh10b_CP000948.1");
	int edit_distance = sv_compute_edit_distance(pure_insert, e_coli_ref);
	ASSERT_EQ(edit_distance, 50); // One insert, one delete, edit distance = 3.5 + 3.5
}

#include <gtest/gtest.h>

#include "modules/bio_format/vcf.h"

TEST(entropy, basic)
{
	std::string seq1 {
		"TGCTCTGAAAAGAGTAACGCGCTTTACTATCCCTGACAATCACCAACAACATCGAACAAGAT"
		"AATAAATTCCTGGTTTAATATCCGACAAGTGAAAACATGCACCCGGACGGGCAGCATGTCGCTCCACAAGTGCAGAGCTT"
		"ACTTGTGTTGTACCGAAGCACTCTGTTCAGGTGGCTGATAGTTGTCAATGTGACTCGCCACGCCAAGAAGAATGACTGAA"
		"ACGACAAGAACGATCCAACCTGTTAATTCAATAAGACGATTCATTACAGCCCACATCTCTCTTGATTGATCCATTAACTT"
		"CAGGGGGTAAATGTTACTTAGCAATAATAGCTCCAGCAAGATTTTTACTGAGGTTTTTGCGATATTTAGCTTTTGTCGTT"
		"GGAAAATTCGCTATTTTTTTGACTTAAGTTAAACAACATCCCTTATTGCTGGCAGGTTATTAAACTGTTGAGCGTGGGTA"
		"AAGGATAGTGTCAAATAGCCATCATACTTCAATGAGAGGCAATGACATGAGCGACAACATCCGTGTTGGGTTGATTGGGTATGG"
	};
	ASSERT_LT(vcf_exporter::compute_entropy(seq1) - 1.978941, 1e-6);
	
	std::string seq2 { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" };
	ASSERT_EQ(vcf_exporter::compute_entropy(seq2), 0.0);

	std::string seq3 { "AAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCC" };
	ASSERT_EQ(vcf_exporter::compute_entropy(seq3), 1.0);
}

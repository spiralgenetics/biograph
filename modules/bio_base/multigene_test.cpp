
#include <gtest/gtest.h>
#include "modules/bio_base/align_multigene.h"
#include "modules/bio_base/align_astar.h"

TEST(multigene, t1)
{
	dna_sequence g1("ACTTACGTAGCTAGCTCAGCTTTAGC");
	dna_sequence g2("CCGTAGAAAACTGACCTGACTAGCTA");
	dna_sequence r("TAGCAGCTCAGAAAACTGACCCTGA");
	dna_sequence r2("TAGCAGCTTAGAAAACTGACCCTGA");
	std::vector<align_info> ai;
	align_multigene(r, g1, g2, ai);
	print_multigene(r, g1, g2, ai);
	align_multigene(r2, g1, g2, ai);
	print_multigene(r2, g1, g2, ai);
}

TEST(multigene, t2)
{
	std::vector<dna_sequence> genes;
	genes.push_back(dna_sequence("ACTTACGTAGCTAGCTCAGCTTTAGC"));
	genes.push_back(dna_sequence("CCGTAGAAAACTGACCTGACTAGCTA"));
	dna_sequence r("TTACTAGCTGACCTGACTAGC");
	std::vector<align_state> path;
	cost_matrix cm;
	cm.ins = 1.5;
	cm.del = 1.5;
	cm.mismatch = 1.0;
	double err = align_astar_skip(path, r, genes, cm, 2.1);
	printf("er = %f\n", err);
	for(size_t i = 0; i < path.size(); i++)
	{
		printf("%d:%d:%d\n", (int) path[i].read_pos, path[i].seq_num, (int) path[i].seq_pos);
	}
}


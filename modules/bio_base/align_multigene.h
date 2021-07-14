
#ifndef __align_multigene_h__
#define __align_multigene_h__

#include "modules/bio_base/dna_sequence.h"

struct align_info
{
	int seq; // Which sequence (0 or 1)
	int pos; // Location in sequence, -1 for deletion
};

double align_multigene(const dna_sequence& r, const dna_sequence& g1, const dna_sequence& g2, std::vector<align_info>& out);
void print_multigene(const dna_sequence& r, const dna_sequence& g1, const dna_sequence& g2, const std::vector<align_info>& out, bool all = true);

#endif

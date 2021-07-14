
#pragma once
#include "modules/bio_base/struct_var.h"

struct allele {
	std::vector<int> sub_ids;
	dna_sequence seq;
	std::vector<size_t> depth;
	std::vector<size_t> fwd;
	std::vector<size_t> tot_qual;
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(sub_ids);
		FIELD(seq);
		FIELD(depth);
		FIELD(fwd);
		FIELD(tot_qual);
	}
};

struct sv_call {
	seq_position position;
	std::vector<struct_var> sources;
	std::vector<allele> alleles;
	double sv_ref_depth; // Ref depth for cis-chromosomal ++ svs
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(position);
		FIELD(sources);
		FIELD(alleles);
		FIELD(sv_ref_depth);
	}
};

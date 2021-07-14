#pragma once

#include <vector>
#include <boost/operators.hpp>
#include "modules/bio_base/align_astar.h"
#include "modules/bio_base/seqset.h"

struct seqset_align_state 
	: boost::less_than_comparable<seqset_align_state>
	, boost::equality_comparable<seqset_align_state>
{
	int read_pos; 
	seqset_range seqset_pos;
	int state;  // 0 = normal, 1 = in ins, 2 = in del, 3 = done
	seqset_align_state(int _read_pos, seqset_range _seqset_pos, int _state = 0) :
		read_pos(_read_pos), seqset_pos(_seqset_pos), state(_state) {}
	bool operator<(const seqset_align_state& rhs) const;
	bool operator==(const seqset_align_state& rhs) const;
};


// Aligns read against sequence, must use all of both
double align_seqset(std::vector<seqset_align_state>& out, const dna_sequence& read, const seqset& the_seqset, const cost_matrix& costs, double max_cost);


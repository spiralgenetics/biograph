
#pragma once

#include "modules/bio_base/dna_sequence.h"
#include <boost/operators.hpp>

/*  
	Variant representation in sv_out 

	Example 1, simple variant:

			 012345678   is_structural = false
	Ref: ACTCATCCA   seq_begin = 3 seq_end = 4
			 ===!I====   left_ref = 2  right_ref = 5
	Seq: ACTT.TCCA
			 0123 4567

	Example 2, stuctural variation
			2222222222
			0123456789
			GACTCATCCA    is_structural = true
			====          seq_begin = 4  seq_end = 6
			GACTTGGCAC    left_ref = 23  right_ref = 89
						====
			CCACAAGCAC
			8888888999
			3456789012
*/

struct sv_out
{
	bool is_structural;  // Bool whether left_ref and right_ref are in same strand
  // Anchor drop happens when you have like:
  // chromosome 1:  abcdxefgh
  // chromosome 2:  ABcdyeFGH
  //
  // assembly: ABcdxeFGH
  //
  // At first glance, this looks like a structural variant from "cdxe"
  // on chromosome 1 to "FGH" on chromosome 2.  However, it could be
  // more likely it's a single change from "y" to "x" on chromsome 2.
  //
  // When this case is detected, we ignore the match to chromosome 1,
  // set the anchor_drop flag, and report it as purely on chromsosome 2.
	bool anchor_drop;    // Is this an anchor drop?

  // If the cost gets too high, we give up on a*, set this flag, and
  // just output a big change.
	bool align_failed;

	int seq_begin;       // Begining of seq var
	int seq_end;         // End of seq var
	dna_const_iterator left_ref;  // Last match on left ref
	dna_const_iterator right_ref; // First match on right ref
	sv_out flip(size_t asm_size) const;
};

struct compound_cost :
	boost::less_than_comparable<compound_cost,
	boost::equality_comparable<compound_cost,
	boost::addable<compound_cost>>>
{
	compound_cost():sv_cost(0), left_alignment_cost(0) {}

	compound_cost(const int _sv_cost, const unsigned long long _left_alignment_cost):
		sv_cost(_sv_cost), left_alignment_cost(_left_alignment_cost)
	{}

	compound_cost(const int _sv_cost):
		sv_cost(_sv_cost), left_alignment_cost(0)
	{}

	void operator+=(const compound_cost& rhs)
	{
		sv_cost += rhs.sv_cost;
		left_alignment_cost += rhs.left_alignment_cost;
	}

	bool operator<(const compound_cost& rhs) const
	{
		if( sv_cost == rhs.sv_cost )
			return left_alignment_cost < rhs.left_alignment_cost;
		return  sv_cost < rhs.sv_cost;
	}

	bool operator==(const compound_cost& rhs) const
	{
		return (sv_cost == rhs.sv_cost) && (left_alignment_cost == rhs.left_alignment_cost);
	}

	int sv_cost;
	unsigned long long left_alignment_cost;
};

struct sv_costs
{
	sv_costs(int _sv_cost = 50, int _var_cost = 4, int _base_cost = 1, int _drop_anchor_cost = 40)
		: sv_cost(_sv_cost), var_cost(_var_cost), base_cost(_base_cost), drop_anchor_cost(_drop_anchor_cost)
	{}
	int sv_cost;
	int var_cost;
	int base_cost;
	int drop_anchor_cost;
};

std::vector<sv_out> call_structural(
	const dna_sequence& seq, 
	const dna_slice& left, 
	const dna_slice& right, 
	compound_cost max_dist,
	const sv_costs& cost = sv_costs()
);



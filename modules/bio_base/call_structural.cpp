
#include "modules/bio_base/call_structural.h"
#include "modules/bio_base/astar.h"
#include "modules/io/log.h"

sv_out sv_out::flip(size_t asm_size) const {
	sv_out out;
	out.is_structural = is_structural;
	out.anchor_drop = anchor_drop;
	out.align_failed = align_failed;
	out.seq_begin = asm_size - seq_end;
	out.seq_end = asm_size - seq_begin;
	out.left_ref = right_ref.rev_comp();
	out.right_ref = left_ref.rev_comp();
	return out;
}

struct sv_astar_loc 
	: boost::less_than_comparable<sv_astar_loc>
	, boost::equality_comparable<sv_astar_loc>
{
	bool in_left;  // Am I currently attched to the left anchor
	bool in_var;    // Am I currently in a variant
	int seq_pos;	// Which is the next base in the seq to process
	dna_const_iterator ref_pos;  // Which is the next ref base to process
	bool operator<(const sv_astar_loc& rhs) const {
		if (in_left != rhs.in_left) return in_left < rhs.in_left;
		if (in_var != rhs.in_var) return in_var < rhs.in_var;
		if (seq_pos != rhs.seq_pos) return seq_pos < rhs.seq_pos;
		if (ref_pos.is_rev_comp() != rhs.ref_pos.is_rev_comp()) return ref_pos.is_rev_comp();
		return ref_pos < rhs.ref_pos;
	}
	bool operator==(const sv_astar_loc& rhs) const {
		return in_left == rhs.in_left &&
			in_var == rhs.in_var && 
			seq_pos == rhs.seq_pos &&
			ref_pos == rhs.ref_pos;
	}
};

/* Starting state

       left_begin     left_end
       |              |
       |              |      seq_end
       v              v      |
       ATACTCGGACTCCGA.......v
       ATACTCAAAGCTATCCGCGCAT
       ^       ACCTATCCGCGCAT
       |       ^             ^
    seq_begin  |             |
               right_begin   right_end

*/

struct sv_astar_context 
{
	compound_cost sv_cost;
	int var_cost;
	int base_cost;
	int drop_anchor_cost;
	const dna_sequence& seq;
	dna_const_iterator left_begin;  
	dna_const_iterator left_end;  
	dna_const_iterator right_begin;
	dna_const_iterator right_end; 
	size_t right_size;  // right_end - right_begin

	typedef compound_cost dist_t;
	typedef sv_astar_loc location_t;

	sv_astar_context(const dna_sequence& _seq, const dna_slice& left, const dna_slice& right, const sv_costs& cost)
		: sv_cost(cost.sv_cost)
		, var_cost(cost.var_cost)
		, base_cost(cost.base_cost)
		, drop_anchor_cost(cost.drop_anchor_cost)
		, seq(_seq)
		, left_begin(left.begin())
		, left_end(left.end())
		, right_begin(right.begin())
		, right_end(right.end())
		, right_size(right.size())
	{}

	compound_cost estimate(const sv_astar_loc& loc, const sv_astar_loc& the_end) const
	{
		//SPLOG("Estimating in_left = %d, in_var = %d, seq_pos = %d, ref_pos = %d", 
		//	loc.in_left, loc.in_var, loc.seq_pos, 
		//	(int) (loc.in_left ? loc.ref_pos - left_begin : loc.ref_pos - right_end));
		int remaining = seq.size() - loc.seq_pos;
		int min_cost = 1000000; 
		// Check best cost to end of right hand side
		if (right_end.is_rev_comp() == loc.ref_pos.is_rev_comp() && // Are we on the same strand
			(loc.ref_pos - right_begin >= 0 && right_end - loc.ref_pos >= 0))  // Are we in range
		{
			int right_remaining = right_end - loc.ref_pos;
			//SPLOG("right_remaining = %d", right_remaining);
			min_cost = base_cost * abs(right_remaining - remaining);
			if (min_cost > 0 && !loc.in_var)
				min_cost += var_cost;
		}
		if (loc.in_left) {
			min_cost = std::min(drop_anchor_cost, min_cost);
		}
		//SPLOG("estimate = %d", min_cost);
		return min_cost;
	}		

	std::vector<std::pair<compound_cost, sv_astar_loc> > nearby(const sv_astar_loc& loc) const
	{
		//SPLOG("Forwarding: in_left = %d, in_var = %d, seq_pos = %d, ref_pos = %d", 
		//	loc.in_left, loc.in_var, loc.seq_pos, 
		//	(int) (loc.in_left ? loc.ref_pos - left_begin : loc.ref_pos - right_end));
		std::vector<std::pair<compound_cost, sv_astar_loc> > out;
		size_t remaining = seq.size() - loc.seq_pos;
		if (loc.in_left && !loc.in_var && remaining <= right_size) { // Can we jump
			sv_astar_loc nl = loc;  // Set up post jump state
			nl.in_left = false;
			dna_const_iterator new_it = right_end - remaining;  // Get *best* spot to jump
			if (remaining == 0) {  
				//SPLOG("Dropping right anchor");
				nl.ref_pos = right_end;
				out.push_back(std::make_pair(drop_anchor_cost, nl));
			} else if (loc.seq_pos == 0) {
				//SPLOG("Dropping left anchor");
				nl.ref_pos = new_it;
				out.push_back(std::make_pair(drop_anchor_cost, nl));
			}
			else if (new_it == loc.ref_pos) { // Is the jump a perfect match
				//SPLOG("Free jump!");
				out.push_back(std::make_pair(0, nl));  // Jump for free!
			} else {
				nl.in_var = true;  // No need to pay twice
				nl.ref_pos = new_it;
				out.push_back(std::make_pair(sv_cost, nl));
			}
		}
		if (loc.in_var && remaining > 0) { // Can we slide seq
			sv_astar_loc nl = loc; 
			nl.seq_pos++;  // Move seq forward
			out.push_back(std::make_pair(base_cost, nl));
		}
		bool is_fwd_safe = (loc.in_left ? (loc.ref_pos != left_end) : (loc.ref_pos != right_end));
		if (loc.in_var && is_fwd_safe) { // Can we slide ref
			sv_astar_loc nl = loc;
			nl.ref_pos++; 
			out.push_back(std::make_pair(base_cost, nl));
		}
		if (loc.in_var && is_fwd_safe && remaining > 0) { // Can we do a mismatch?
			sv_astar_loc nl = loc;
			nl.seq_pos++;
			nl.ref_pos++;
			out.push_back(std::make_pair(base_cost, nl));
		}
		if (!loc.in_var && is_fwd_safe && remaining > 0 && 
			seq[loc.seq_pos] == *loc.ref_pos) { // Can we do a match?
			sv_astar_loc nl = loc;
			nl.seq_pos++;
			nl.ref_pos++;
			out.push_back(std::make_pair(0, nl));  // Free!
		}
		if (loc.in_var) { // Try ending var
			sv_astar_loc nl = loc;
			nl.in_var = false;
			out.push_back(std::make_pair(0, nl));
		}
		if (!loc.in_var) { // Try starting var
			sv_astar_loc nl = loc;
			nl.in_var = true;
			// note here that we add the left-alignment bias:
			// the farther the reference pos for the start of the var, the bigger the penalty.
			out.push_back(std::make_pair(compound_cost(var_cost, loc.ref_pos.get_original_offset()), nl));
		}
		if (!loc.in_left && loc.in_var && loc.ref_pos != right_begin) { // Special left moving case
			sv_astar_loc nl = loc;
			nl.ref_pos--;
			out.push_back(std::make_pair(base_cost, nl));
		}
		return out;
	}
};

std::vector<sv_out> call_structural(
	const dna_sequence& seq, 
	const dna_slice& left, 
	const dna_slice& right, 
	compound_cost max_dist, 
	const sv_costs& cost)
{
	//SPLOG("Doing structural call, seq = %s", seq.as_string().c_str());
	std::vector<sv_out> out;
	sv_astar_context context(seq, left, right, cost);
	sv_astar_loc start;
	sv_astar_loc end;
	start.in_left = true;
	start.in_var = false;
	start.seq_pos = 0;
	start.ref_pos = left.begin();
	end.in_left = false;
	end.in_var = false;
	end.seq_pos = seq.size();
	end.ref_pos = right.end();
	astar_state<sv_astar_context> as(context, start, end, max_dist);
	compound_cost the_cost = as.run();
	if (the_cost == max_dist)
	{
		//SPLOG("Gave up due to high cost, doing simple version");
		sv_out simple;
		simple.is_structural = true;
		simple.align_failed = true;
		simple.seq_begin = 0;
		simple.left_ref = left.begin();
		simple.seq_end = seq.size();
		simple.right_ref = right.end() - 1;
		while( simple.left_ref != left.end() && 
			simple.seq_begin < (int) seq.size() &&
			*simple.left_ref == seq[simple.seq_begin])
		{
			simple.left_ref++;
			simple.seq_begin++;
		}
		simple.left_ref--;
		while( simple.right_ref + 1 != right.begin() &&
			simple.seq_end > 0 &&
			simple.seq_begin != simple.seq_end &&
			*simple.right_ref == seq[simple.seq_end - 1])
		{
			simple.right_ref --;
			simple.seq_end--;
		}
		simple.right_ref++;
		//if (simple.seq_begin != simple.seq_end) 
		out.push_back(simple);
		return out;
	}
	//SPLOG("Final cost = %d", the_cost)
	std::vector<sv_astar_loc> path;
	as.get_path(path);
	//SPLOG("Generating path")
	
	sv_out cur;
	cur.align_failed = false;
	for(size_t i = 1; i < path.size(); i++)
	{
		sv_astar_loc prev = path[i-1];
		sv_astar_loc next = path[i];
		//SPLOG("in_left = %d, in_var = %d, seq_pos = %d (%c), ref_pos = %d (%c)", 
		//	next.in_left, next.in_var, next.seq_pos, (char) seq[next.seq_pos],
		//	(int) (next.in_left ? next.ref_pos - left.begin : next.ref_pos - right.end),
		//	(char) *next.ref_pos);
		if (prev.in_left == true && next.in_left == false)
		{
			if (prev.seq_pos == (int) seq.size())
			{
				//printf("RIGHT ANCHOR DROP: %s\n", seq.as_string().c_str());
				cur.is_structural = false;
				cur.anchor_drop = true;
				cur.seq_begin = cur.seq_end = seq.size();
				cur.left_ref = prev.ref_pos;
				cur.right_ref = next.ref_pos;
				out.push_back(cur);
			}
			else if (prev.seq_pos == 0) {
				//printf("LEFT ANCHOR DROP: %s\n", seq.as_string().c_str());
				cur.is_structural = false;
				cur.anchor_drop = true;
				cur.seq_begin = cur.seq_end = 0;
				cur.left_ref = prev.ref_pos;
				cur.right_ref = next.ref_pos;
				out.push_back(cur);
			}
			else if (prev.ref_pos != next.ref_pos) 
			{
				//SPLOG("Structural Variation START!");
				cur.is_structural = true;
				cur.anchor_drop = false;
				cur.seq_begin = next.seq_pos;
				cur.left_ref = prev.ref_pos - 1;
			}
		} else if (!prev.in_var && next.in_var) {
			//SPLOG("Normal Variation START!");
			cur.is_structural = false;
			cur.anchor_drop = false;
			cur.seq_begin = next.seq_pos;
			cur.left_ref = next.ref_pos - 1;
		} else if (prev.in_var && !next.in_var) {
			//SPLOG("Variation END!");
			cur.seq_end = next.seq_pos;
			cur.right_ref = next.ref_pos;
			out.push_back(cur);
		}
	}

	return out;
}



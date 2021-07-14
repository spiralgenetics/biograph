#include "modules/bio_base/astar.h"
#include "modules/bio_mapred/align_seqset.h"

bool seqset_align_state::operator<(const seqset_align_state& rhs) const
{
	if (read_pos != rhs.read_pos)
		return read_pos < rhs.read_pos;
	if (state != rhs.state)
		return state < rhs.state;
	return seqset_pos < rhs.seqset_pos;
}

bool seqset_align_state::operator==(const seqset_align_state& rhs) const
{
	return read_pos == rhs.read_pos && 
		seqset_pos == rhs.seqset_pos &&
		state == rhs.state;
}

struct seqset_align_astar_state
{
	typedef seqset_align_state location_t;
	typedef double dist_t;

	const seqset& the_seqset;
	const dna_sequence& read;
	const cost_matrix& costs;

	seqset_align_astar_state(const seqset& _the_pwbt, const dna_sequence& _read, const cost_matrix& _costs)
		: the_seqset(_the_pwbt)
		, read(_read)
		, costs(_costs)	
	{}

	double estimate(const seqset_align_state& a, const seqset_align_state& b) const { return 0.0; }

	// Finds nearby states
	std::vector<std::pair<double, seqset_align_state>> nearby(const seqset_align_state& loc) const
	{
		std::vector<std::pair<double, seqset_align_state>> out;
		if (loc.state == 3) { 
			// If we are already done, nowhere to go
			return out; 
		}
		if (loc.read_pos == (int) read.size()) {
			// If we are at the end of read, move to done
			out.emplace_back(0.0, seqset_align_state(read.size(), the_seqset.end(), 3));
			return out;
		}

		// Try adding each base
		for(int b = 0; b < 4; b++) {
			seqset_range next = loc.seqset_pos.push_front(dna_base(b));
			if (!next.valid()) {
				// If there is nowhere to go, skip
				continue;
			}
			if (next.size() > 70) {
				next = next.pop_back();
			}
			// Handle non-insert case
			out.emplace_back(
				// Zero costs if it's a match
				(b == (int) read[loc.read_pos].complement()) ? 0.0 : costs.mismatch,
				seqset_align_state(loc.read_pos + 1, next, 0));
			// Handle insert case
			out.emplace_back(
				loc.state == 1 ? costs.extend_ins : costs.ins,
				seqset_align_state(loc.read_pos, next, 1));
		}
		// Handle delete case
		out.emplace_back(loc.state == 2 ? costs.extend_del : costs.del,
				seqset_align_state(loc.read_pos + 1, loc.seqset_pos, 2));
		return out;
	}
};

double align_seqset(std::vector<seqset_align_state>& out, const dna_sequence& read, const seqset& the_seqset, const cost_matrix& costs, double max_cost)
{
	seqset_align_astar_state ctx(the_seqset, read, costs);
	seqset_align_state start(0, the_seqset.ctx_begin(), 0 );
	seqset_align_state goal(read.size(), the_seqset.end(), 3);
	astar_state<seqset_align_astar_state> as(ctx, start, goal, max_cost);
	double r = as.run();
	if (r < max_cost)
		as.get_path(out);
	return r;
}



#include "modules/bio_base/align_astar.h"
#include "modules/bio_base/astar.h"

bool align_state::operator<(const align_state& rhs) const
{
	if (read_pos != rhs.read_pos)
		return read_pos < rhs.read_pos;
	if (seq_num != rhs.seq_num)
		return seq_num < rhs.seq_num;
	return seq_pos < rhs.seq_pos;
}

bool align_state::operator==(const align_state& rhs) const
{
	return read_pos == rhs.read_pos && 
		seq_num == rhs.seq_num &&
		seq_pos == rhs.seq_pos;
}

struct align_astar_base
{
	typedef align_state location_t;
	typedef double dist_t;

	align_astar_base(const dna_sequence& _read, const cost_matrix& _costs)
		: read(_read)
		, costs(_costs)	
	{}

	const dna_sequence& read;
	const cost_matrix& costs;

	double estimate(const align_state& a, const align_state& b) const { return 0.0; }

	// Computes simple changes (not jumping from seq to seq)
	void simple_nearby(std::vector<std::pair<double, align_state> >& out, const dna_sequence& seq, const align_state& loc) const
	{
		// Get current sequence I'm on
		if (loc.read_pos < (int) read.size() && loc.seq_pos < (int) seq.size())  // If both seqs can more forward
		{	
			// Try match case (with zero cost)
			if (seq[loc.seq_pos] == read[loc.read_pos])
				out.push_back(std::make_pair(0, 
					align_state(loc.read_pos + 1, loc.seq_num, loc.seq_pos + 1)));
			else  // Mismatch case
				out.push_back(std::make_pair(costs.mismatch, 
					align_state(loc.read_pos + 1, loc.seq_num, loc.seq_pos + 1)));
		}
		
		if (loc.read_pos < (int) read.size()) // If not on end of read
		{
			// Handle insert case
			out.push_back(std::make_pair((loc.last_op == 1 ? costs.extend_ins : costs.ins),
				align_state(loc.read_pos + 1, loc.seq_num, loc.seq_pos, 1)));
		}
		if (loc.seq_pos < (int) seq.size())  // If not at end of seq
		{
			// Handle delete case
			out.push_back(std::make_pair((loc.last_op == 2 ? costs.extend_del : costs.del),
				align_state(loc.read_pos, loc.seq_num, loc.seq_pos + 1, 2)));
		}
	}
};

template<class ctx_t>
double wrap_as_func(std::vector<align_state>& out, const dna_sequence& read, typename ctx_t::seq_t seq, const cost_matrix& costs, double max_cost)
{
	ctx_t ctx(read, seq, costs);
	align_state start = ctx.get_start();
	align_state goal = ctx.get_goal(); 
	astar_state<ctx_t> as(ctx, start, goal, max_cost);
	double r = as.run();
	if (r < max_cost)
		as.get_path(out);
	return r;
}

struct align_astar_exact_context : public align_astar_base
{
	typedef const dna_sequence& seq_t;
	const dna_sequence& seq;
	align_astar_exact_context(const dna_sequence& _read, const dna_sequence& _seq, const cost_matrix& _costs)
		: align_astar_base(_read, _costs), seq(_seq) {}

	std::vector<std::pair<double, align_state> > nearby(const align_state& loc) const
	{
		std::vector<std::pair<double, align_state> > out;
		simple_nearby(out, seq, loc);
		return out;
	}
	align_state get_start() { return align_state(0, 0, 0); }
	align_state get_goal() { return align_state(read.size(), 0, seq.size()); }
};

double align_astar_exact(std::vector<align_state>& out, const dna_sequence& read, const dna_sequence& seq, const cost_matrix& costs, double max_cost)
{
	return wrap_as_func<align_astar_exact_context>(out, read, seq, costs, max_cost);
}

struct align_astar_float_context : public align_astar_base
{
	typedef const dna_sequence& seq_t;
	const dna_sequence& seq;
	align_astar_float_context(const dna_sequence& _read, const dna_sequence& _seq, const cost_matrix& _costs)
		: align_astar_base(_read, _costs), seq(_seq) {}

	std::vector<std::pair<double, align_state> > nearby(const align_state& loc) const
	{
		std::vector<std::pair<double, align_state> > out;
		if (loc.seq_num == 1) 
			return out;
		if (loc.seq_pos == -1)
		{
			align_state s = loc;
			for(size_t i = 0; i < seq.size(); i++)
			{
				s.seq_pos = i;
				out.push_back(std::make_pair(0, s));
			}
		}
		else 
		{
			simple_nearby(out, seq, loc);
			if (loc.read_pos == (int) read.size())
				out.push_back(std::make_pair(0,
					align_state(loc.read_pos, 1, -1)));
		}
		return out;
	}
	align_state get_start() { return align_state(0, 0, -1); }
	align_state get_goal() { return align_state(read.size(), 1, -1); }
};

double align_astar_float(std::vector<align_state>& out, const dna_sequence& read, const dna_sequence& seq, const cost_matrix& costs, double max_cost)
{
	return wrap_as_func<align_astar_float_context>(out, read, seq, costs, max_cost);
}
struct align_astar_skip_context : public align_astar_base
{
	typedef const std::vector<dna_sequence>& seq_t;
	seq_t seqs;
	align_astar_skip_context(const dna_sequence& _read, seq_t _seqs, const cost_matrix& _costs)
		: align_astar_base(_read, _costs), seqs(_seqs) {}

	std::vector<std::pair<double, align_state> > nearby(const align_state& loc) const
	{
		std::vector<std::pair<double, align_state> > out;
		if (loc.seq_num == (int) seqs.size()) 
			return out;
		if (loc.seq_pos == -1)
		{
			align_state s = loc;
			for(size_t i = 0; i < seqs[loc.seq_num].size(); i++)
			{
				s.seq_pos = i;
				out.push_back(std::make_pair(0, s));
			}
		}
		else 
		{
			simple_nearby(out, seqs[loc.seq_num], loc);
			out.push_back(std::make_pair(0,
				align_state(loc.read_pos, loc.seq_num + 1, -1)));
		}
		return out;
	}
	align_state get_start() { return align_state(0, 0, -1); }
	align_state get_goal() { return align_state(read.size(), seqs.size(), -1); }
};

double align_astar_skip(std::vector<align_state>& out, const dna_sequence& read, const std::vector<dna_sequence>& seqs, const cost_matrix& costs, double max_cost)
{
	return wrap_as_func<align_astar_skip_context>(out, read, seqs, costs, max_cost);
}

struct align_astar_any_context : public align_astar_base
{
	typedef const std::vector<dna_sequence>& seq_t;
	seq_t seqs;
	align_astar_any_context(const dna_sequence& _read, seq_t _seqs, const cost_matrix& _costs)
		: align_astar_base(_read, _costs), seqs(_seqs) {}

	std::vector<std::pair<double, align_state> > nearby(const align_state& loc) const
	{
		std::vector<std::pair<double, align_state> > out;
		if (loc.seq_num == (int) seqs.size()) 
			return out;
		if (loc.seq_pos == -1)
		{
			for(size_t i = 0; i < seqs.size(); i++)
			{
				for(size_t j = 0; j < seqs[i].size(); j++)
				{
					out.push_back(std::make_pair(0, 
						align_state(loc.read_pos, i, j)));
				}
			}
		}
		else 
		{
			simple_nearby(out, seqs[loc.seq_num], loc);
			if (loc.read_pos == (int) read.size())
				out.push_back(std::make_pair(0,
					align_state(loc.read_pos, seqs.size(), -1)));
		}
		return out;
	}
	align_state get_start() { return align_state(0, 0, -1); }
	align_state get_goal() { return align_state(read.size(), seqs.size(), -1); }
};

double align_astar_any(std::vector<align_state>& out, const dna_sequence& read, const std::vector<dna_sequence>& seqs, const cost_matrix& costs, double max_cost)
{
	return wrap_as_func<align_astar_any_context>(out, read, seqs, costs, max_cost);
}



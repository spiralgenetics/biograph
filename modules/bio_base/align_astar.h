
#ifndef __align_astar_h__
#define __align_astar_h__

#include <vector>
#include <boost/operators.hpp>
#include "modules/bio_base/dna_sequence.h"

struct align_state 
	: boost::less_than_comparable<align_state>
	, boost::equality_comparable<align_state>
{
	int read_pos; 
	int seq_num; 
	int seq_pos;  // Pos of binding, -1 if still unattached
	int last_op;  // 0 = normal, 1 = ins, 2 = del
	align_state(int _read_pos, int _seq_num, int _seq_pos, int _last_op = 0) :
		read_pos(_read_pos), seq_num(_seq_num), seq_pos(_seq_pos), last_op(_last_op) {}
	bool operator<(const align_state& rhs) const;
	bool operator==(const align_state& rhs) const;
};

namespace std
{
	template<>
	inline void swap(align_state& a, align_state& b) noexcept
	{
		swap(a.read_pos, b.read_pos);
		swap(a.seq_num, b.seq_num);
		swap(a.seq_pos, b.seq_pos);
	}
}

struct cost_matrix
{
	cost_matrix() : ins(3.5), del(3.5), mismatch(1.0), extend_ins(1.5), extend_del(1.5) {}
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(ins, TF_STRICT);
		FIELD(del, TF_STRICT);
		FIELD(mismatch, TF_STRICT);
		FIELD(extend_ins, 1.5);
		FIELD(extend_del, 1.5);
	}
	double ins;
	double del;
	double mismatch;
	double extend_ins;
	double extend_del;
};

// Aligns read against sequence, must use all of both
double align_astar_exact(std::vector<align_state>& out, const dna_sequence& read, const dna_sequence& seq, const cost_matrix& costs, double max_cost);

// Aligns read against sequence, all of read must be accounted for, but not seq
double align_astar_float(std::vector<align_state>& out, const dna_sequence& read, const dna_sequence& seq, const cost_matrix& costs, double max_cost);

// Aligns read against sequences in order, skipping as optimal
double align_astar_skip(std::vector<align_state>& out, const dna_sequence& read, const std::vector<dna_sequence>& seqs, const cost_matrix& costs, double max_cost);

// Aligns read agains any of seqs (or none) as optimal 
double align_astar_any(std::vector<align_state>& out, const dna_sequence& read, const std::vector<dna_sequence>& seqs, const cost_matrix& costs, double max_cost);

#endif


#include "base/base.h"
#include "modules/bio_base/pileup.h"
#include <queue>
#include "modules/io/log.h"

pileup::pileup(const dna_sequence& sequence, size_t max_cost)
	: m_sequence(sequence)
	, m_max_cost(max_cost)
{
	m_pileup.resize(m_sequence.size());
}

int pileup::add_read(const std::string& name, const dna_sequence& read_seq,
	 const std::string quality, bool fwd, int offset)
{
	// Find where the read belongs if needed
	if (offset == -1) {
		offset = best_match_pos(read_seq, quality);
		if (offset < 0) return offset;
	}
	// Add to pileup
	for(size_t i = 0; i < read_seq.size(); i++)
	{
		base_info& bi = m_pileup[i + offset][read_seq[i]];
		bi.count++;
		bi.tot_qual += uint32_t(quality[i]);
		if (fwd) bi.fwd++;
	}

	return offset;
}

size_t pileup::depth_at(size_t position) const
{
	if (position >= m_sequence.size())
	{
		SPLOG("BAD DEPTH: Attempting to get depth at position %lu", position);
		return 0;
	}
	auto it = m_pileup[position].find(m_sequence[position]);
	if (it == m_pileup[position].end()) { return 0; }
	return it->second.count;
}

size_t pileup::fwd_at(size_t position) const
{
	auto it = m_pileup[position].find(m_sequence[position]);
	if (it == m_pileup[position].end()) { return 0; }
	return it->second.fwd;
}

size_t pileup::tot_qual_at(size_t position) const
{
	auto it = m_pileup[position].find(m_sequence[position]);
	if (it == m_pileup[position].end()) { return 0; }
	return it->second.tot_qual;
}

/*
std::vector<pileup::read_info> pileup::relevant_reads(size_t start, size_t end)
{
	std::vector<read_info> out;
	return out;
}
*/

// Returns position of best match or -1 if no match
int pileup::best_match_pos(const dna_sequence& read_seq, const std::string quality)
{
	// Setup priority queue
	struct pq_entry 
	{ 
		size_t score; int seq_pos; int read_pos; 
		pq_entry(size_t _score, int _seq_pos, int _read_pos) 
			: score(_score), seq_pos(_seq_pos), read_pos(_read_pos) {}
		bool operator<(const pq_entry& rhs) const { return score > rhs.score; }
	};
	std::priority_queue<pq_entry> queue;

	// If seq is too long, forget it
	if (read_seq.size() > m_sequence.size()) return -1;

	// Add initial entries
	for(size_t i = 0; i <= m_sequence.size() - read_seq.size(); i++) 
	{
		queue.push(pq_entry(0, i, 0));
	}

	while(!queue.empty()) {
		// Pop top entry
		pq_entry cur = queue.top();
		queue.pop();
		// If too large a cost, give up
		if (cur.score > m_max_cost) return -1;
		// If fully matched, return result
		if (size_t(cur.read_pos) == read_seq.size())
			return cur.seq_pos - cur.read_pos;
		// Add cost if mismatch
		if (read_seq[cur.read_pos] != m_sequence[cur.seq_pos])
			cur.score += size_t(quality[cur.read_pos]);
		// Move positions forward
		cur.read_pos++;
		cur.seq_pos++;
		queue.push(cur);
	}
	// Should never set here
	CHECK(false);
	return -1;
}

void pileup::print()
{
	for(size_t i = 0; i < m_sequence.size(); i++)
	{
		std::string bases;
		for(int bn = 0; bn < 4; bn++)
		{
			dna_base b(bn);
			int count = m_pileup[i][b].count;
			if (count > 0)
				bases += printstring("%c(%d) ", (char) b, count);
		}
		SPLOG("%d: %s", (int) i, bases.c_str());
	}
}


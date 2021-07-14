
#include "modules/bio_base/align_multigene.h"
#include "modules/io/io.h"
#include <boost/multi_array.hpp>

static double match_score(const dna_base& b1, const dna_base& b2)
{
	if (b1 == b2) return 0.0;
	return 1.0;
}

static const double insert_score = 2.5;  // Read has 'insert' not in original
static const double delete_score = 2.5;  // Read skips over base in original

struct choice
{
	choice() : score(0.0), which(-1) {}
	choice(double _score) : score(_score), which(-1) {}
	choice(double _score, int _which) : score(_score), which(_which) {}
	double score;
	int which;
};

choice choose(const choice& c1, const choice& c2)
{
	if (c2.score < c1.score)
		return c2;
	return c1;
}

choice choose(const choice& c1, const choice& c2, const choice& c3)
{ return choose(c1, choose(c2, c3)); }

choice choose(const choice& c1, const choice& c2, const choice& c3, const choice& c4)
{ return choose(c1, choose(c2, c3, c4)); }

double align_multigene(const dna_sequence& r, const dna_sequence& g1, const dna_sequence& g2, std::vector<align_info>& out)
{
	boost::multi_array<choice, 2> rg1(boost::extents[r.size()+1][g1.size()]);
	boost::multi_array<choice, 2> rg2(boost::extents[r.size()+1][g2.size()]);
	boost::multi_array<choice, 1> rs(boost::extents[r.size()+1]); 
	// It's free to start 'r' anywhere
	for(size_t i = 0; i < g1.size(); i++)
		rg1[0][i] = 0;
	// Score for skipping parts of r
	for(size_t i = 1; i <= r.size(); i++)
		rg1[i][0] = insert_score * i;
	// Main body of rg1
	for(size_t i = 1; i <= r.size(); i++)
	{
		for(size_t j = 1; j < g1.size(); j++)
		{
			rg1[i][j] = choose(
				choice(rg1[i-1][j-1].score + match_score(r[i-1], g1[j-1]), 1),
				choice(rg1[i][j-1].score + delete_score, 2),
				choice(rg1[i-1][j].score + insert_score, 3)
				);
		}
	}
	// Calculate cost of various skip points
	for(size_t i = 0; i <= r.size(); i++)
	{
		rs[i] = choice(rg1[i][0].score, 0);
		for(size_t j = 1; j < g1.size(); j++)
			rs[i] = choose(rs[i], choice(rg1[i][j].score, j));
	}
	// Score of first alignment anywhere in g2
	for(size_t i = 0; i < g2.size(); i++)
		rg2[0][i] = 0;
	// Score carried across from other side 
	for(size_t i = 1; i <= r.size(); i++)
		rg2[i][0] = choice(rs[i].score, 0);
	// Main body of rg2
	for(size_t i = 1; i <= r.size(); i++)
	{
		for(size_t j = 1; j < g2.size(); j++)
		{
			rg2[i][j] = choose(
				choice(rs[i].score, 0),
				choice(rg2[i-1][j-1].score + match_score(r[i-1], g2[j-1]), 1),
				choice(rg2[i][j-1].score + delete_score, 2),
				choice(rg2[i-1][j].score + insert_score, 3)
				);
		}
	}
	// Find the minimal result that fully consumes the read
	choice final(rg2[r.size()][0].score, 0);
	for(size_t i = 1; i < g2.size(); i++)
		final = choose(final, choice(rg2[r.size()][i].score, i));

	// Now travel backwards and generate output
	out.resize(r.size());
	int cur_seq = 1;
	int cur_pos = final.which;
	int cur_r = r.size();
	while(cur_r)
	{
		choice back = (cur_seq ? rg2[cur_r][cur_pos] : rg1[cur_r][cur_pos]);
		//printf("cur_seq = %d, cur_pos = %d, cur_r = %d, back.score = %f, back.which = %d\n",
		//	(int) cur_seq, (int) cur_pos, (int) cur_r, back.score, (int) back.which);
		switch(back.which)
		{
		case 0:
			cur_seq = 0;
			cur_pos = rs[cur_r].which;
			break;
		case 1:
			cur_r--;
			cur_pos--;
			out[cur_r].seq = cur_seq;
			out[cur_r].pos = cur_pos;
			break;
		case 2:
			cur_pos--;
			break;
		case 3:
			cur_r--;
			out[cur_r].seq = cur_seq;
			out[cur_r].pos = -1;
			break;
		default:
			throw io_exception("Invalid case in back.which");  // Should never happen, but just in case
		}
	}
	return final.score;
}

void print_multigene(const dna_sequence& r, const dna_sequence& g1, const dna_sequence& g2, const std::vector<align_info>& out, bool all)
{
	int cur_seq = 0;
	int cur_pos = 0;
	std::string topline;
	std::string botline;
	bool started = false;
	for(size_t i = 0; i < r.size(); i++)
	{
		// Deal with deletions
		if (out[i].pos == -1)
		{
			started = true;
			topline += '.';
			botline += char(r[i]);
			continue;
		}
		// Print sequence 1 bases missing in read
		while(cur_seq == 0 && (out[i].seq == 1 || (out[i].seq == 0 && cur_pos < out[i].pos)))
		{
			if (all || started)
			{
				topline += char(g1[cur_pos]);
				botline += ' ';
				if (!all)
				{
					topline += '-';
					botline += ' ';
					started = false;
				}
			}
			cur_pos++;
			if (cur_pos == (int) g1.size()) {cur_seq = 1; cur_pos = 0; topline += ' '; botline += ' ';}
		}
		// Print sequence 2 bases missing in read
		while(cur_seq == 1 && out[i].seq == 1 && cur_pos < out[i].pos)
		{
			if (all || started)
			{
				topline += char(g2[cur_pos]);
				botline += ' ';
				if (!all)
				{
					topline += '-';
					botline += ' ';
					started = false;
				}
			}
			cur_pos++;
		}
		// Print 'matched' base
		topline += char((out[i].seq ? g2[cur_pos] : g1[cur_pos]));
		cur_pos++;
		botline += char(r[i]);
		started = true;
	}
	if (all)
	{
		while(cur_seq == 0)
		{
			topline += char(g1[cur_pos++]);
			if (cur_pos == (int) g1.size()) {cur_seq = 1; cur_pos = 0; topline += ' ';}
		}
		while(cur_seq == 1 && cur_pos != (int) g2.size())
			topline += char(g2[cur_pos++]);
	}
	printf("%s\n%s\n", topline.c_str(), botline.c_str());
}


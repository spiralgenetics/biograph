
#include "modules/bio_base/dna_multiseq.h"
#include <boost/multi_array.hpp>

// 4 diffs is better than 1 insert/delete, but 5 diffs is worse
const int del_cost = 9;
const int diff_cost = 4;

dna_multiseq::dna_multiseq(const dna_sequence& s1, const dna_sequence& s2)
{
	int s1s = s1.size();
	int s2s = s2.size();
	boost::multi_array<int, 2> cost(boost::extents[s1.size()+1][s2.size()+1]);
	boost::multi_array<int, 2> dir(boost::extents[s1.size()+1][s2.size()+1]);
	cost[s1s][s2s] = 0;
	for(int i = 0; i < s1s; i++)
	{
		cost[i][s2s] = del_cost;
		dir[i][s2s] = 1;
	}
	for(int j = 1; j < s2s; j++)
	{
		cost[s1s][j] = del_cost;
		dir[s1s][j] = 2;
	}
	for(int i = s1s-1; i >= 0; i--)
	{
		for(int j = s2s-1; j >= 0; j--)
		{
			int mcost = cost[i+1][j+1] + (s1[i] != s2[j] ? diff_cost : 0);
			int d1cost = cost[i+1][j] + del_cost;
			int d2cost = cost[i][j+1] + del_cost;
			if (mcost <= d1cost && mcost <= d2cost)
			{
				dir[i][j] = 0;
				cost[i][j] = mcost;
			}
			else if (d1cost <= d2cost)
			{
				dir[i][j] = 1;
				cost[i][j] = d1cost;
			}
			else
			{
				dir[i][j] = 2;
				cost[i][j] = d2cost;
			}
		}
	}	

	int i = 0;
	int j = 0;
	
	m_seqs.resize(2);
	dna_del_seq& o1 = m_seqs[0];
	dna_del_seq& o2 = m_seqs[1];
	
	while(i < s1s || j < s2s)
	{
		if (dir[i][j] == 0)
		{
			o1.push_back(dna_del_base((char) s1[i]));
			o2.push_back(dna_del_base((char) s2[j]));
			i++; j++;
		}
		else if (dir[i][j] == 1)
		{
			o1.push_back(dna_del_base((char) s1[i]));
			o2.push_back(dna_del_base('.'));
			i++;
		}
		else 
		{
			o1.push_back(dna_del_base('.'));
			o2.push_back(dna_del_base((char) s2[j]));
			j++;
		}
	}
}

std::string dna_multiseq::get_string(size_t which)
{
	dna_del_seq& seq = m_seqs[which];
	std::string out;
	for(size_t i = 0; i < seq.size(); i++)
		out += (char) seq[i];
	return out;
}


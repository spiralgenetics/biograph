#include "modules/bio_base/kmer.h"
#include "modules/bio_mapred/align_kmer.h"
#include "modules/bio_base/astar.h"
#include "modules/io/log.h"
#include <boost/operators.hpp>

static int best_pow(int loc)
{
	for(int i = 0; i < 16; i++)
	{
		if ((loc & (1 << i)) != 0)
			return i;
	}
	return 16;
}

struct kmer_astar_state
	: boost::less_than_comparable<kmer_astar_state>
	, boost::equality_comparable<kmer_astar_state>
{
	kmer_astar_state(int _start, int _end, kmer_t _ks, kmer_t _ke)
		: start(_start)
		, end(_end)
		, ks(_ks)
		, ke(_ke)
	{
		pow2s = best_pow(start);
		pow2e = best_pow(end);
		min_pow = std::min(pow2s, pow2e);
		max_pow = std::max(pow2s, pow2e);
	}
		
	int start;
	int end;
	kmer_t ks;
	kmer_t ke;
	int pow2s;
	int pow2e;
	int min_pow;
	int max_pow;
	bool operator<(const kmer_astar_state& rhs) const
	{
		if (min_pow != rhs.min_pow)
			return min_pow < rhs.min_pow;
		if (max_pow != rhs.max_pow)
			return max_pow < rhs.max_pow;
		if (start != rhs.start)
			return start < rhs.start;
		if (end != rhs.end)
			return end < rhs.end;
		return false;
	}
	bool operator==(const kmer_astar_state& rhs) const
	{
		return (start == rhs.start && end == rhs.end);
	}
};

struct kmer_astar_context
{
	typedef kmer_astar_state location_t;
	typedef double dist_t;
	kmer_astar_context(const dna_sequence& read, const std::string& qual, const kmer_set& kmers, double min_base_quality)
		: m_read(read)
		, m_qual(qual)
		, m_kmers(kmers)
		, m_min_base_quality(min_base_quality)
		, m_ks(kmers.kmer_size())
	{}

	double estimate(const kmer_astar_state& a, const kmer_astar_state& b) const { return 0.0; }

	std::vector<std::pair<double, kmer_astar_state> > nearby(const kmer_astar_state& loc) const
	{
		std::vector<std::pair<double, kmer_astar_state> > r;
		if (loc.start == loc.end) 
		{
			// Initial state, find all matching kmers
			for(size_t i = 0; i <= m_read.size() - m_ks; i++)
			{
				kmer_t k = make_kmer(m_read.begin() + i, m_ks);
				if (m_kmers.count(canonicalize(k, m_ks)))
					r.push_back(std::make_pair(0.0,
						kmer_astar_state(i, i + m_ks, k, k)));
			}
		}
		else
		{
			//SPLOG("Processing rs=%d (%d:%d) %s:%s", m_read.size(), loc.start, loc.end, 
			//	dna_sequence(loc.ks, m_ks).as_string().c_str(),
			//	dna_sequence(loc.ke, m_ks).as_string().c_str());
			if (loc.pow2s < loc.pow2e && loc.start > 0)
			{
				// Try adding in reverse direction
				kmer_t k = loc.ks;
				for(int base = 0; base < 4; base++)
				{
					kmer_t k2 = append(base, left(k, m_ks, m_ks-1), m_ks-1);
					if (m_kmers.count(canonicalize(k2, m_ks)))
					{
						double cost = (dna_base(base) == m_read[loc.start - 1]) ? 0.0 
							: std::max(m_min_base_quality, double(m_qual[loc.start - 1] - 33));
						r.push_back(std::make_pair(cost,
							kmer_astar_state(loc.start - 1, loc.end, k2, loc.ke)));
					}
				}
			}
			else if (loc.end < (int) m_read.size())
			{
				// Try adding in forward direction
				kmer_t k = loc.ke;
				for(int base = 0; base < 4; base++)
				{
					kmer_t k2 = append(right(k, m_ks-1), base, 1);
					if (m_kmers.count(canonicalize(k2, m_ks)))
					{
						double cost = (dna_base(base) == m_read[loc.end]) ? 0.0 
							: std::max(m_min_base_quality, double(m_qual[loc.end] - 33));
						r.push_back(std::make_pair(cost,
							kmer_astar_state(loc.start, loc.end +1, loc.ks, k2)));
					}
				}
			}
			else if (loc.start == 0 && loc.end == (int) m_read.size() && loc.ks && loc.ke)
			{
				r.push_back(std::make_pair(0.0, kmer_astar_state(0, m_read.size() + 1, 0, 0)));
			}
		}
		return r;
	}

	const dna_sequence& m_read;
	std::string m_qual;
	const kmer_set& m_kmers;
	double m_min_base_quality;
	size_t m_ks;
};

unsigned verify_kmers(const dna_sequence& read, const kmer_set& kmers) {
  size_t kmer_size = kmers.kmer_size();
  if (read.size() < kmer_size) return 0;
  for (size_t i = 0; i <= read.size() - kmer_size; i++) {
    kmer_t k = make_kmer(read.begin() + i, kmer_size);
    k = canonicalize(k, kmer_size);
    if (!kmers.count(k)) {
      return i + kmer_size - 1;
    }
  }
  return read.size();
  return true;
}

double align_kmer(std::vector<kmer_t>& out, const dna_sequence& read, const std::string& qual, 
		const kmer_set& kmers, double min_base_quality, double max_cost)
{
	if (read.size() < kmers.kmer_size()) return max_cost;
	
	kmer_astar_context ctx(read, qual, kmers, min_base_quality);
	kmer_astar_state start(0, 0, 0, 0);
	kmer_astar_state end(0, read.size() + 1, 0, 0);
	astar_state<kmer_astar_context> as(ctx, start, end, max_cost);
	double cost = as.run();
	if (cost >= max_cost)
		return max_cost;
	std::vector<kmer_astar_state> r;
	as.get_path(r);
	out.resize(read.size() - kmers.kmer_size() + 1);
	for(size_t i = 0; i < r.size(); i++)
	{
		if (r[i].end == 0 || r[i].end > (int) read.size()) continue;
		//printf("Adding (%d, %d) = (%s, %s)\n", r[i].start, r[i].end, 
		//		dna_sequence(r[i].ks, kmers.kmer_size()).as_string().c_str(), dna_sequence(r[i].ke, kmers.kmer_size()).as_string().c_str());
		out[r[i].start] = r[i].ks;
		out[r[i].end - kmers.kmer_size()] = r[i].ke;
	}
	return cost;
}

dna_sequence get_corrected(const std::vector<kmer_t>& in, size_t kmer_size)
{
	dna_sequence start(in[0], kmer_size);
	for(size_t i = 1; i < in.size(); i++)
	{
		start.push_back(dna_base((int) right(in[i], 1)));
	}
	return start;
}


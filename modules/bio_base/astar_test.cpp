
#include <map>
#include <gtest/gtest.h>
#include "modules/bio_base/astar.h"

int best_pow(int loc)
{
	for(int i = 0; i < 16; i++)
	{
		if ((loc & (1 << i)) != 0)
			return i;
	}
	return 16;
}

struct mypair
{
	mypair(int _start, int _end)
		: start(_start)
		, end(_end)
	{
		pow2s = best_pow(start);
		pow2e = best_pow(end);
		min_pow = std::min(pow2s, pow2e);
		max_pow = std::max(pow2s, pow2e);
	}
		
	int start;
	int end;
	int pow2s;
	int pow2e;
	int min_pow;
	int max_pow;
	bool operator<(const mypair& rhs) const
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
	bool operator==(const mypair& rhs) const
	{
		return (start == rhs.start && end == rhs.end);
	}
};

mypair expand(const mypair& mp)
{
	if (mp.pow2s < mp.pow2e && mp.start > 0)
		return mypair(mp.start - 1, mp.end);
	else
		return mypair(mp.start, mp.end + 1);
}

struct test_astar_context
{
	typedef mypair location_t;
	typedef double dist_t;
	double estimate(const mypair& a, const mypair& b) const { return 0.0; }
	std::vector<std::pair<double, mypair> > nearby(const mypair& loc) const
	{
		std::vector<std::pair<double, mypair> >  r;
		if (loc.start == loc.end) 
		{
			for(int i = 3; i < 17; i++)
				r.push_back(std::make_pair(0.0, mypair(i, i + 7)));
		}
		else
		{
			mypair out = expand(loc);
			if (out.end <= 27)
			{
				printf("extending (%d, %d) -> (%d, %d)\n", loc.start, loc.end, out.start, out.end);
				r.push_back(std::make_pair(0.0, expand(loc)));
			}
		}
		return r;
	}
};

TEST(astar, test1)
{
	test_astar_context ctx;
	mypair start(0, 0);
	mypair end(0, 27);
	astar_state<test_astar_context> as(ctx, start, end, 10.0);
	as.run();
}


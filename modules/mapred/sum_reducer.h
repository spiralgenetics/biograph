
#ifndef __sum_reducer_h__

#include "modules/mapred/reducer.h"
#include "modules/io/double_pair.h"


class sum_reducer : public simple_reducer<sum_reducer, uint64_t>
{
public:
	sum_reducer(const std::string& params) {}

	virtual void typed_summarize(uint64_t& total, const uint64_t& add) { total += add; }
};


class pair_sum_reducer : public simple_reducer<pair_sum_reducer, double_pair>
{
public:
	pair_sum_reducer(const std::string& params) {}

	virtual void typed_summarize(double_pair& total, const double_pair& add) { total += add; }
};


#endif


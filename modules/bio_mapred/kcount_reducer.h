#pragma once

#include "modules/mapred/reducer.h"
#include "modules/bio_base/kmer.h"

class kcount_reducer : public simple_reducer<kcount_reducer, kcount_pair>
{
public:
	kcount_reducer(const std::string& params) {}

	virtual void typed_summarize(kcount_pair& total, const kcount_pair& add) {
		total.fwd += add.fwd;
		total.rev += add.rev;
	}
};


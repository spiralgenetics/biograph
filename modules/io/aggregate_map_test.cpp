
#define AGGREGATE_MAP_DEBUG

#include "base/base.h"
#include "gtest/gtest.h"
#include "modules/io/aggregate_map.h"
#include <map>

const size_t k_op_count = 5000;
const size_t k_key_range = 100000;
const size_t k_value_range = 100;

TEST(AGGREGATE_MAP, basic)
{
	printf("Testing aggregate_map\n");
	aggregate_map<int, int> t;
	std::map<int, int> dt;

	for(size_t i = 0; i < k_op_count; i++)
	{
		if (i % 1000 == 0)
			printf("Op %d\n", (int) i);
		CHECK(t.size() == dt.size());
		switch(random() % 5)
		{
		case 0:
		case 1:
		{
			int k = random() % k_key_range;
			aggregate_map<int, int>::const_iterator it = t.find(k);
			std::map<int, int>::const_iterator it2 = dt.find(k);
			if (it2 == dt.end())
			{
				CHECK(it == t.end());
				int v = random() % k_value_range;
				//printf("Doing insert of %d\n", k);
				t.insert(std::make_pair(k, v));
				dt.insert(std::make_pair(k, v));
			}
			else
				CHECK_EQ(it->second, it2->second);
		}
		break;
		case 2:
		{
			int k = random() % k_key_range;
			std::map<int, int>::const_iterator it = dt.upper_bound(k);
			if (it != dt.end())
			{
				t.erase(it->first);
				dt.erase(it->first);
			}
		}
		break;
		case 3:
		{
			int k1 = random() % k_key_range;
			int k2 = random() % k_key_range;
			if (k2 < k1)
				std::swap(k1, k2);
			//printf("Checking range (%d-%d)\n", k1, k2);
			std::map<int, int>::const_iterator it = dt.lower_bound(k1);
			std::map<int, int>::const_iterator itEnd = dt.lower_bound(k2);
			int total = 0;
			for(; it != itEnd; ++it)
				total += it->second;
			int logtot = t.total(k1, k2);
			//printf("total = %d, logtot = %d\n", total, logtot);
			if (total != logtot)
			{
				t.dump();
				exit(1);
			}
			CHECK_EQ(total, logtot);
		}
		case 4:
		{
			int k = random() % k_key_range;
			aggregate_map<int, int>::const_iterator it = t.lower_bound(k);
			std::map<int, int>::const_iterator it2 = dt.lower_bound(k);
			if (it == t.end())
				CHECK(it2 == dt.end());
			else
				CHECK_EQ(it->first, it2->first);
		}	
		break;
		}
		//t.dump();
		std::map<int, int>::const_iterator it;
		aggregate_map<int, int>::const_iterator it2 = t.begin();
		for(it = dt.begin(); it != dt.end(); ++it, ++it2)
		{
			CHECK_EQ(it->first, it2->first);
			CHECK_EQ(it->second, it2->second);
		}
		CHECK(it2 == t.end());
			
		CHECK(t.validate());
	}
}

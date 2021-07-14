
#include "modules/io/prefix_sum.h"
#include "modules/io/log.h"
#include <vector>
#include <gtest/gtest.h>
#include <memory>

class fake_prefix_sum
{
public:
	fake_prefix_sum(size_t size) : m_vec(size, 0) {}
	void add(size_t which, size_t val) { m_vec[which] += val; }
	size_t total(size_t which) { 
		size_t tot = 0;
		for(size_t i = 0; i < which; i++) {
			tot += m_vec[i];
		}
		return tot;
	}		
	void nearest_below(uint32_t x, uint32_t& idx, uint32_t& tot) {
		tot = 0;
		idx = 0;
		while(tot + m_vec[idx] <= x) {
			tot += m_vec[idx];
			idx++;
		}
	}
private:
	std::vector<size_t> m_vec;
};

TEST(prefix_sum, test)
{
	size_t ps_size = 121;
	prefix_sum ps(ps_size);
	fake_prefix_sum fps(ps_size);
	// Put in one initial value so things aren't ill-formed
	ps.add(50, 1);
	fps.add(50, 1);
	for(size_t i = 0; i < 100000; i++) {
		int op = random() % 3;
		if (op == 0) {
			// Add to some value
			size_t which = random()%ps_size;
			size_t val = random()%5;
			fps.add(which, val);
			ps.add(which, val);
		} else if (op == 1) {
			// Check totals
			size_t which = random() % (ps_size + 1);
			ASSERT_EQ(fps.total(which), ps.total(which));
		} else if (op == 2) {
			uint32_t num = random() % ps.total();
			uint32_t idx, tot;
			uint32_t fidx, ftot;
			ps.nearest_below(num, idx, tot);
			fps.nearest_below(num, fidx, ftot);
			ASSERT_LE(ps.total(idx), num);
			ASSERT_GT(ps.total(idx + 1), num);
			ASSERT_EQ(idx, fidx);
			ASSERT_EQ(tot, ftot);
		}
	}
}

TEST(prefix_sum, grow)
{
	prefix_sum the_prefix_sum;
	
	ASSERT_EQ(the_prefix_sum.size(), 0);
	
	for (unsigned i = 0; i < 1024; i++) {
		the_prefix_sum.push_back(i);
		ASSERT_EQ(the_prefix_sum.size(), i + 1);
		ASSERT_EQ(the_prefix_sum.total(0), 0);
		ASSERT_EQ(the_prefix_sum.total(the_prefix_sum.size() - 1), i * (i - 1) / 2);
}	}
	

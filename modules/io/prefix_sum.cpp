#include "modules/io/log.h"

#include "modules/io/prefix_sum.h"

static uint32_t next_greatest_power_of_2(uint32_t in) 
{
	uint32_t pot = 1;
	while(pot < in) { pot *= 2; }
	return pot;
}

prefix_sum::prefix_sum(uint32_t size)
	: m_size(size)
	, m_potsize(next_greatest_power_of_2(size))
	, m_tree(2*m_potsize, 0)
{}

void prefix_sum::add(uint32_t which, uint32_t val)
{
	uint32_t cur = m_potsize + which;
	while(cur) {
		m_tree[cur] += val;
		cur /= 2;
	}
}

void prefix_sum::sub(uint32_t which, uint32_t val)
{
	uint32_t cur = m_potsize + which;
	while(cur) {
		m_tree[cur] -= val;
		cur /= 2;
	}
}

uint32_t prefix_sum::total(uint32_t which) const
{
	if (which == m_potsize) return m_tree[1];
	uint32_t cur = m_potsize + which;
	uint32_t tot = 0;
	while(cur != 1) {
		tot += (cur & 1) * m_tree[cur-1];
		cur /= 2;
	}
	return tot;
}

void prefix_sum::nearest_below(uint32_t x, uint32_t& idx, uint32_t& tot) const
{
	uint32_t cur = 1;
	tot = 0;
	while(cur < m_potsize) {
		if (tot + m_tree[2*cur] <= x) {
			tot += m_tree[2*cur];
			cur = 2*cur + 1;
		} else {
			cur = 2*cur;
		}
	}
	idx = cur - m_potsize;
}

void prefix_sum::push_back(uint32_t new_value)
{
	if (m_size == m_potsize) {
		uint32_t original_potsize = m_potsize;
		m_potsize *= 2;
		m_tree.resize(2 * m_potsize);
		while (original_potsize) {
			std::swap_ranges(m_tree.begin() + original_potsize
				, m_tree.begin() + original_potsize * 2
				, m_tree.begin() + original_potsize * 2
			);
			original_potsize /= 2;
		}
		m_tree[1] = m_tree[2];
	}
	
	add(m_size++, new_value);
}

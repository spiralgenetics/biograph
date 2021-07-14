#pragma once

#include "modules/io/io.h"

class prefix_sum
{
public:
	explicit prefix_sum(uint32_t size);
	prefix_sum() : prefix_sum(0) {}

	void reset() { std::fill(m_tree.begin(), m_tree.end(), 0U); }

	void add(uint32_t which, uint32_t val);
	void sub(uint32_t which, uint32_t val);
	uint32_t total(uint32_t which) const;
	uint32_t total() const { return m_tree[1]; }
	
	void push_back(uint32_t new_value);
	uint32_t size() const { return m_size; }
	
	uint32_t value(uint32_t which) const { return m_tree[m_potsize + which]; }
	// Find the last index such that total(idx) <= x, returns idx, tot= total(idx)
	void nearest_below(uint32_t x, uint32_t& idx, uint32_t& tot) const;

private:
	uint32_t m_size;
	uint32_t m_potsize; // Power of two size.
	std::vector<uint32_t> m_tree;
};

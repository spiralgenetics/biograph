
#pragma once

#include <vector>
#include <tuple>
#include <stdint.h>
#include <math.h>

#include "base/base.h"
#include "modules/io/io.h"
#include "modules/io/log.h"
#include "modules/io/prefix_sum.h"

typedef uint32_t symbol_t;
typedef uint32_t range_t;
typedef uint32_t count_t;

constexpr range_t k_max_range = 0x80000000;
constexpr range_t k_half_range = (k_max_range / 2);
constexpr range_t k_min_range = (k_half_range / 2);

/*
0xxx:0yyy  -> 0, xxx0:yyy0
1xxx:1yyy  -> 1, xxx0:yyy0
1xxx:0yyy  -> Impossible
00xx:11yy  -> Good!
01xx:10yy  -> 0xx0:1yy0
*/
inline int rc_increase_range(range_t& start, range_t& end)
{
	CHECK_LT(start, end);
	if (start >= k_half_range)
	{
		start = (start - k_half_range) * 2;
		end = (end - k_half_range) * 2;
		return 1;
	}
	else if (end < k_half_range)
	{
		start *= 2;
		end *= 2;
		return -1;
	}
	else
	{
		// 01xx:10yy  -> 0xx0:1yy0
		// CHECK((start & k_half_range) == 0);
		// CHECK((start & k_min_range) == k_min_range);
		// CHECK((end & k_half_range) == k_half_range);
		// CHECK((end & k_min_range) == 0);
		start = (start - k_min_range) * 2;
		end = (end - k_min_range) * 2;
		return 0;
	}
}

class range_encoder
{
public:
	range_encoder(writable& out)
		: m_out(out)
		, m_start(0)
		, m_end(k_max_range)
		, m_deferred(0)
		, m_byte(0)
		, m_bit_count(0)
		, m_entropy(0.0)
		, m_entropy_stats(false)
	{}

	template<class Model>
	void encode(const Model& model, symbol_t s)
	{
		// Subtend the range based on the symbol
		range_t len = m_end - m_start;

		//SPLOG("old range = [%x, %x), symbol = %d, ", m_start, m_end, s);
		uint32_t start, end;
		model.symbol_range(s, len, start, end);
		m_end = m_start + end;
		m_start += start;
		//SPLOG("new range = [%x, %x), symbol = %d, ", m_start, m_end, s);

		// This is expensive, so only collect stats if they were requested.
		if(m_entropy_stats) {
			m_entropy += -log2(double(m_end - m_start) / double(len));
		}

		while(m_end - m_start < k_min_range)
		{
			switch(rc_increase_range(m_start, m_end))
			{
			case -1:
				//SPLOG("0: start = %x, end = %x\n", m_start, m_end);
				put(0);
				for(count_t i = 0; i < m_deferred; i++)
					put(1);
				m_deferred = 0;
				break;
			case 0:
				//SPLOG("D: start = %x, end = %x\n", m_start, m_end);
				m_deferred++;
				break;
			case 1:
				//SPLOG("1: start = %x, end = %x\n", m_start, m_end);
				put(1);
				for(count_t i = 0; i < m_deferred; i++)
					put(0);
				m_deferred = 0;
				break;
			default:
				CHECK(false);
			}
		}
	}

	void end()
	{
		while(m_start != 0)
		{
			if (m_start >= k_half_range)
			{
				put(1);
				for(count_t i = 0; i < m_deferred; i++)
					put(0);
				m_deferred = 0;
				m_start = (m_start - k_half_range) * 2;
			}
			else
			{
				put(0);
				for(count_t i = 0; i < m_deferred; i++)
					put(1);
				m_deferred = 0;
				m_start *= 2;
			}
		}
		if (m_deferred)
		{
			put(1);
		}
		if (m_bit_count == 0)
			return;
		m_byte <<= (8 - m_bit_count);
		m_out.write((const char*) &m_byte, 1);
	}

	double cur_entropy() { return m_entropy; }
	void reset_entropy() { m_entropy = 0.0; }
	void enable_stats() { m_entropy_stats = true; }

private:
	void put(symbol_t bit)
	{
		CHECK_LT(bit, 2);
		m_byte <<= 1;
		m_byte |= bit;
		m_bit_count++;
		if (m_bit_count == 8)
		{
			m_out.write((const char*) &m_byte, 1);
			m_byte = 0;
			m_bit_count = 0;
		}
	}
	writable& m_out;
	range_t m_start;
	range_t m_end;
	count_t m_deferred;
	uint8_t m_byte;
	uint8_t m_bit_count;
	double  m_entropy;
	bool m_entropy_stats;
};

class range_decoder
{
public:
	range_decoder(readable& in)
		: m_in(in)
		, m_start(0)
		, m_end(1)
		, m_val(0)
		, m_byte(0)
		, m_bit_count(0)
	{
		while(m_end < k_max_range)
		{
			m_end <<= 1;
			m_val <<= 1;
			m_val |= (range_t) get();
		}
	}

	template<class Model>
	symbol_t decode(const Model& model)
	{
		// Regenerate the symbol from a binary search
		range_t len = m_end - m_start;
		symbol_t s;
		range_t start;
		range_t end;
		//SPLOG("old range = [%x, %x)", m_start, m_end);
		model.symbol_find(m_val - m_start, len, s, start, end);
		m_end = m_start + end;
		m_start += start;
		//SPLOG("new range = [%x, %x), symbol = %d, ", m_start, m_end, s);

		while(m_end - m_start < k_min_range)
		{
			if (rc_increase_range(m_start, m_end) == 0)
				m_val -= k_min_range;
			if (m_val >= k_half_range)
				m_val -= k_half_range;
			m_val <<= 1;
			m_val |= (range_t) get();
		}

		return s;
	}

private:
	template<class Model>
	inline symbol_t binary_search(Model& model, range_t cur, range_t len)
	{
		symbol_t bot = 0;
		symbol_t top = model.count();

		while(top - bot != 1)
		{
			symbol_t mid = bot + (top - bot) / 2;
			if (model.range_start(mid, len) <= cur) {
				bot = mid;
			} else {
				top = mid;
			}
		}
		return bot;
	}

	symbol_t get()
	{
		if (m_bit_count == 0)
		{
			m_bit_count = 8;
			m_byte = 0;
			// Ignore EOF and generate 0's
			m_in.read((char*) &m_byte, 1);
		}
		m_bit_count--;
		return (m_byte & (1 << m_bit_count) ? 1 : 0);
	}

private:
	readable& m_in;
	range_t m_start;
	range_t m_end;
	range_t m_val;
	uint8_t m_byte;
	uint8_t m_bit_count;
};

class uniform_dist
{
public:
	uniform_dist(symbol_t count) : m_count(count) {}
	symbol_t count() const { return m_count; }
	void symbol_range(symbol_t s, range_t r, range_t& start, range_t& end) const
	{
		start = uint64_t(r) * uint64_t(s) / uint64_t(m_count);
		end = uint64_t(r) * uint64_t(s + 1) / uint64_t(m_count);
	}
	void symbol_find(range_t x, range_t r, symbol_t& s, range_t& start, range_t& end) const
	{
		s = (uint64_t(x + 1) * uint64_t(m_count) - uint64_t(1)) / uint64_t(r);
		symbol_range(s, r, start, end);
	}
private:
	symbol_t m_count;
};

class prefix_sum_dist
{
public:
	prefix_sum_dist() : m_ps() {}
	prefix_sum_dist(size_t size) : m_ps(size) {}
	symbol_t count() const { return m_ps.size(); }
	void symbol_range(symbol_t s, range_t r, range_t& start, range_t& end) const
	{
		uint64_t grand_tot = m_ps.total();
		uint64_t start_tot = m_ps.total(s);
		uint64_t end_tot = start_tot + m_ps.value(s);
		start = start_tot * uint64_t(r) / grand_tot;
		end = end_tot * uint64_t(r) / grand_tot;
		//SPLOG("Symbol = %u, Range = %u, orig = [%d,%d), out=[%u,%u)", s, r,
		//	(int) start_tot, (int) end_tot, start, end);
	}
	void symbol_find(range_t x, range_t r, symbol_t& s, range_t& start, range_t& end) const
	{
		uint64_t grand_tot = m_ps.total();
		uint32_t seek_off = (uint64_t(x + 1) * uint64_t(grand_tot) - uint64_t(1)) / uint64_t(r);
		uint32_t tot;
		m_ps.nearest_below(seek_off, s, tot);
		uint64_t start_tot = tot;
		uint64_t end_tot = start_tot + m_ps.value(s);
		start = start_tot * uint64_t(r) / grand_tot;
		end = end_tot * uint64_t(r) / grand_tot;
		//SPLOG("x=%u, r = %u, seek = %u, Symbol = %u, out=[%u,%u)",
		//	x, r, seek_off, s, start, end);
	}
	prefix_sum& inner() { return m_ps; }
	const prefix_sum& inner() const { return m_ps; }
private:
	prefix_sum m_ps;
};



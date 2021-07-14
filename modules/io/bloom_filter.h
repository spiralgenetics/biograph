#pragma once

#include "modules/io/packed_vector.h"

#include <cmath>

template <size_t CounterWidth>
class bloom_filter
{
public:
	typedef mutable_packed_vector<size_t, CounterWidth> bitmap_type;

	struct traits
	{
		traits() 
			: cells(0)
			, hashes(0)
		{}

		explicit traits(size_t capacity, double error_rate)
			: cells(compute_cells(capacity, error_rate))
			, hashes(compute_hashes(capacity, cells))
		{}

		size_t cells;
		size_t hashes;

	private:
		static const size_t compute_hashes(size_t capacity, size_t cells)
		{
			return ceil(cells * log(2) / capacity);
		}

		static const size_t compute_cells(size_t capacity, double error_rate)
		{
			return ceil(capacity * fabs(log(error_rate)) / pow(log(2), 2));
		}
	};

	explicit bloom_filter(const traits& traits)
		: m_traits(traits)
		, m_bitmap(m_traits.cells, "bloom_filter")
	{}

	explicit bloom_filter(size_t capacity, double error_rate)
		: bloom_filter(traits(capacity, error_rate))
	{}

public:
	size_t cells() const 
	{
		return m_traits.cells;
	}

	size_t hashes() const
	{
		return m_traits.hashes;
	}

	const bitmap_type& bitmap() const
	{
      return m_bitmap;
	}

public:
	void clear()
	{
		m_bitmap.reset();
	}

	template <typename Hasher>
	bool add(const Hasher& hasher)
	{
		bool overflow = false;
		for (size_t i = 0; i < m_traits.hashes; i++) {
			auto digest = hasher(i);
			auto index = digest % m_bitmap.size();
			overflow |= m_bitmap[index].safe_increment();
		}
		return overflow;
	}

	template <typename Hasher>
	size_t lookup(const Hasher& hasher) const
	{
		size_t count = m_bitmap.max_value()+1;
		for (size_t i = 0; i < m_traits.hashes; i++) {
			auto digest = hasher(i);
			auto index = digest % m_bitmap.size();
			auto cell = m_bitmap[index];
			count = std::min(count, (size_t)cell);
		}
		return count;
	}

private:
	traits m_traits;
	bitmap_type m_bitmap;
};

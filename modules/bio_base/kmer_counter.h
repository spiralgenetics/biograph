#pragma once

#include "modules/bio_base/kmer.h"
#include <vector>
#include <functional>
#include <iterator>
#include <thread>

template <class Hash = std::hash<kmer_t>>
class kmer_counter
{
public:
	struct value_type
	{
		kmer_t key;
		uint32_t fwd_count;
		uint32_t rev_count;
	};

private:
	class iterator : public std::iterator<std::forward_iterator_tag, value_type>
	{
	public:
		explicit iterator(kmer_counter& kc, size_t index)
			: m_kc(kc)
			, m_index(index)
		{
			if (m_index == 0) {
				advance();
			}
		}

		value_type operator *() 
		{
			return value_type{m_kc.m_keys[m_index], m_kc.m_fwd_values[m_index], m_kc.m_rev_values[m_index]};
		}

		iterator& operator ++()
		{
			m_index++;
			advance();
			return *this;
		}
	
		bool operator ==(const iterator& rhs) const { return m_index == rhs.m_index; }
		bool operator !=(const iterator& rhs) const { return m_index != rhs.m_index; }

	private:
		void advance()
		{
			for (; m_index < m_kc.m_keys.size(); m_index++) {
				if (m_kc.m_keys[m_index] == m_kc.k_sentinel) {
					continue;
				}
				break;
			}
		}

	private:
		const kmer_counter& m_kc;
		size_t m_index;
	};

	friend class iterator;

public:
	explicit kmer_counter(size_t capacity)
		: m_keys(capacity, k_sentinel)
		, m_fwd_values(capacity)
		, m_rev_values(capacity)
	{}

public:
	void clear()
	{
		std::fill(m_keys.begin(), m_keys.end(), k_sentinel);
		std::fill(m_fwd_values.begin(), m_fwd_values.end(), 0);
		std::fill(m_rev_values.begin(), m_rev_values.end(), 0);
	}

	size_t capacity() const
	{
		return m_keys.size();
	}

	void add(const kmer_t& key, bool fwd)
	{
		auto index = get_index(key);
		auto probe = (fwd ? &m_fwd_values[index] : &m_rev_values[index]);
		while (true) {
			auto old_value = *probe;
			auto new_value = old_value + 1;
			if (__sync_bool_compare_and_swap(probe, old_value, new_value)) {
				return;
			}
			std::this_thread::yield();
		}
	}

	iterator begin()
	{
		return iterator(*this, 0);
	}

	iterator begin() const
	{
		return iterator(*this, 0);
	}

	iterator end()
	{
		return iterator(*this, m_keys.size());
	}

	iterator end() const
	{
		return iterator(*this, m_keys.size());
	}

private:
	size_t get_index(const kmer_t& key)
	{
		Hash hasher;
		size_t count = 0;
		for (size_t index = hasher(key); count < m_keys.size(); count++) {
			index %= m_keys.size();
			auto probe = &m_keys[index];
			if (*probe != key) {
				if (*probe != k_sentinel) {
					index++;
					continue;
				}

				if (!__sync_bool_compare_and_swap(probe, k_sentinel, key)) {
					continue;
				}
			}
			return index;
		}
		throw io_exception(
			boost::format("The hash table is full, please increase the hash table size, "
				"which is currently set to %lu") % m_fwd_values.capacity()
		);
	}

private:
	const kmer_t k_sentinel = std::numeric_limits<kmer_t>::max();
	std::vector<kmer_t> m_keys;
	std::vector<uint32_t> m_fwd_values;
	std::vector<uint32_t> m_rev_values;
};

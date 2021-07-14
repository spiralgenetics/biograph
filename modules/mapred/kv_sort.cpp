
#include "modules/mapred/kv_sort.h"
#include <algorithm>
#include <boost/foreach.hpp>

kv_sort::kv_pair::kv_pair(const std::string& _key, const std::string& _value)
	: key(_key)
	, value(_value) 
{}

kv_sort::ptr_compare::ptr_compare(std::shared_ptr<sorter> order) 
	: m_order(order) 
{}

bool kv_sort::ptr_compare::operator()(pair_ptr a, pair_ptr b)
{
	int r = m_order->compare(a->key, b->key);
	if (r == 0)
		return a < b;
	return r < 0;
}

kv_sort::kv_sort(const std::string& serialized_params) 
	: m_sorted(true)
	, m_records_size(0)
	, m_order(nullptr)
	, m_current(m_pairs.begin())
{
	kv_sort::param sort_params;
	msgpack_deserialize(sort_params, serialized_params);
	
	m_sorter = sorter_registry::get_safe(sort_params.sorter, "");
	if (! sort_params.splitter.empty())
	{
		m_splitter = splitter_registry::get_safe(sort_params.splitter, sort_params.first_key);
	}
	
	m_order = ptr_compare(m_sorter);
}

kv_sort::~kv_sort()
{
	clear();
}

void kv_sort::write(const std::string& key, const std::string& value)
{
	m_sorted = false;
	pair_ptr pair = new kv_pair(key, value);
	m_pairs.push_back(pair);
	m_records_size += kv_serial_size(key.size(), value.size());
}

void kv_sort::prep_read()
{
	if (!m_sorted)
	{
		std::sort(m_pairs.begin(), m_pairs.end(), m_order);
		m_current = m_pairs.begin();
		m_sorted = true;
	}
}

bool kv_sort::read(std::string& key, std::string& value)
{
	if (m_current == m_pairs.end())
		return false;

	key = (*m_current)->key;
	value = (*m_current)->value;
	m_current++;
	return true;
}

void kv_sort::reset()
{
	m_current = m_pairs.begin();
}

void kv_sort::clear()
{
	BOOST_FOREACH(pair_ptr p, m_pairs)
		delete p;
	m_pairs.clear();
	m_current = m_pairs.begin();
	m_records_size = 0;
}

void kv_sort::set_file_info(file_info& fi) const
{
	if (!m_sorted)
		throw io_exception("kv_sort wasn't sorted before trying to get file_info");
	fi.size = get_size();
	fi.num_records = get_num_records();
	fi.first_key = (*m_pairs.begin())->key;
	fi.last_key = (*m_pairs.rbegin())->key;
}



#include "modules/mapred/kv_summarize.h"
#include <algorithm>
#include <boost/foreach.hpp>

kv_summarize::key_compare::key_compare(const sorter& order) 
	: m_order(order) 
{}

bool kv_summarize::key_compare::operator()(const std::string& a, const std::string& b)
{
	int r = m_order.compare(a, b);
	if (abs(r) < 2)
		return false;  // Things are 'equal'
	return r < 0;
}

kv_summarize::kv_summarize(const std::string& in) 
	: m_params(in)
	, m_records_size(0)
	, m_sorter(sorter_registry::get_safe(m_params.sort, ""))
	, m_reducer(reducer_registry::get_safe(m_params.reduce, m_params.reduce_param))
	, m_order(*m_sorter)
	, m_pairs(m_order)
	, m_current(m_pairs.begin())
{}

kv_summarize::~kv_summarize()
{
	clear();
}

void kv_summarize::write(const std::string& key, const std::string& value)
{
	collection_t::iterator it = m_pairs.find(key);
	if (it == m_pairs.end())
	{
		m_records_size += kv_serial_size(key.size(), value.size());
		m_pairs[key] = value;
	}
	else
	{
		m_records_size -= kv_serial_size(key.size(), it->second.size());
		m_reducer->summarize(it->second, value);
		m_records_size += kv_serial_size(key.size(), it->second.size());
	}
}

void kv_summarize::prep_read()
{
	m_current = m_pairs.begin();
}

bool kv_summarize::read(std::string& key, std::string& value)
{
	if (m_current == m_pairs.end())
		return false;

	key = (*m_current).first;
	value = (*m_current).second;
	m_current++;
	return true;
}

void kv_summarize::reset()
{
	m_current = m_pairs.begin();
}

void kv_summarize::clear()
{
	m_pairs.clear();
	m_current = m_pairs.begin();
	m_records_size = 0;
}

void kv_summarize::set_file_info(file_info& fi) const
{
	fi.size = get_size();
	fi.num_records = get_num_records();
	fi.first_key = m_pairs.begin()->first;
	fi.last_key = m_pairs.rbegin()->first;
}


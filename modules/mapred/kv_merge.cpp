#include "modules/mapred/kv_merge.h"
#include "modules/io/log.h"
#include "modules/mapred/file_info_reader.h"

kv_merge::kv_pair::kv_pair(size_t _part, const std::string& _key, const std::string& _value)
	: part(_part)
	, key(_key)
	, value(_value) 
{}

kv_merge::ptr_compare::ptr_compare(sorter& order) 
	: m_order(order) 
{}

bool kv_merge::ptr_compare::operator()(pair_ptr a, pair_ptr b)
{
	return m_order.compare(a->key, b->key) > 0;
}

kv_merge::kv_merge(sorter& order) 
	: m_queue(ptr_compare(order))
{}

kv_merge::~kv_merge()
{
	while(!m_queue.empty())
	{
		delete m_queue.top();
		m_queue.pop();
	}
}
		
void kv_merge::add(file_info_reader* fir)
{
	size_t part = m_sources.size();
	std::string value;

	m_sources.push_back(fir);
	std::string key = fir->get_first_key();
	if(!key.empty())
	{
		m_queue.push(new kv_pair(part, key, ""));
	}
	else { throw io_exception("trying to merge an unsorted file or a sorted file with no first key set."); }
}

bool kv_merge::read(std::string& key, std::string& value)
{
	hack:
	//SPLOG("Reading from merged file");
	if (m_queue.empty())
		return false;

	pair_ptr pair = m_queue.top();
	m_queue.pop();

	size_t part = pair->part;
	//SPLOG("Found key %s in part %d", pair->key.c_str(), part);

	//is this the very first key?
	if(pair->value == "")
	{
		if (!m_sources[part]->read(key, value)) 
			goto hack;
	}
	else
	{
		key = pair->key;
		value = pair->value;
	}

	if (m_sources[part]->read(pair->key, pair->value))
	{//re-use the pair, instead of doing a new-delete for each read.
		m_queue.push(pair);
	}
	else
	{//only delete the pair when the source is empty.
		delete pair;
	}
	return true;
}


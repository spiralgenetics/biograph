
#include "modules/mapred/merge_reader.h"
#include "modules/io/log.h"

merge_reader::~merge_reader()
{
	delete m_merged;
	for(size_t i = 0; i < m_readers.size(); i++)
		delete m_readers[i];
}

bool merge_reader::read(std::string& key, std::string& value)
{
	//SPLOG("Doing merge read, begin_on = '%s', end_before = '%s', clean_break = %d", 
	//		m_begin_on.c_str(), m_end_before.c_str(), m_clean_break);
	while(m_merged->read(key, value))
	{
		// Only once we are at m_begin_on (or it's group if clean_break) do we start
		if (m_begin_on != "" && !(m_clean_break ?
				m_sorter->compare(m_begin_on, key) < 2 :
				m_sorter->compare(m_begin_on, key) < 1
				)) continue;
		// If we get to m_end_bore (or it's group) we stop
		if (m_end_before != "" && (m_clean_break ? 
				m_sorter->compare(m_end_before, key) < 2 :
				m_sorter->compare(m_end_before, key) < 1
				)) return false;
		return true;
	}

	return false;
}


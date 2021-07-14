
#include "modules/mapred/kv_hold.h"

void kv_hold::write(const std::string& key, const std::string& value)
{
	m_records_size += kv_serial_size(key.size(), value.size());
	m_keys.push_back(key);
	m_values.push_back(value);
}

bool kv_hold::read(std::string& key, std::string& value)
{
	if (m_keys.size() == m_offset)
		return false;

	key = m_keys[m_offset];
	value = m_values[m_offset];
	m_offset++;
	return true;
}

bool kv_hold::legal_split(const std::string& key) const
{
	if (m_keys.size() == 0) return true;
	if (m_sorter == NULL) return true;
	return (abs(m_sorter->compare(key, *m_keys.rbegin())) == 2);
}

void kv_hold::set_file_info(file_info& fi) const 
{ 
	fi.size = get_size();
	// fi.actual_size will be set by whoever does the encoding to the file_info
	fi.num_records = get_num_records(); 
	if (m_sort != "")
	{
		fi.first_key = *m_keys.begin();
		fi.last_key = *m_keys.rbegin();
	}
}

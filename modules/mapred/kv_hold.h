
#ifndef __kv_hold_h__
#define __kv_hold_h__

#include "modules/io/keyvalue.h"
#include "modules/mapred/sorter.h"
#include "modules/mapred/manifest.h"
#include <vector>
#include <queue>

class kv_hold : public kv_sink, public reset_kv_source
{
public:
	kv_hold(const std::string& sort) 
		: m_offset(0), m_records_size(0), m_sort(sort), m_sorter(nullptr),
          m_keys(track_alloc("kv_hold:keys")), m_values(track_alloc("kv_hold:values"))
	{ if (m_sort != "") m_sorter = sorter_registry::get(sort, ""); }
	void write(const std::string& key, const std::string& value) override;
	void prep_read() {} //potentially modify the underlying data before read ops start
	bool read(std::string& key, std::string& value) override;
	void reset() override { m_offset = 0; }
	void clear() { m_keys.clear(); m_values.clear(); m_records_size = 0; m_offset = 0; }
	bool oversized(size_t goal_size) const { return (get_size() + 20 * get_num_records()) > goal_size; } 
	bool legal_split(const std::string& key) const;
	size_t get_num_records() const { return m_keys.size(); }
	size_t get_size() const { return m_records_size; }
	void set_file_info(file_info& fi) const;
	bool split_now(const std::string& key) const { return false; }
	void update_split(const std::string& key) {}

private:
	size_t m_offset;
	size_t m_records_size;
	std::string m_sort;
	std::unique_ptr<sorter>  m_sorter;
	tracked_vector<std::string> m_keys;
	tracked_vector<std::string> m_values;
};

#endif

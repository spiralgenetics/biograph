
#ifndef __merge_reader_h__
#define __merge_reader_h__

#include "modules/mapred/kv_merge.h"
#include "modules/mapred/file_info_reader.h"

// Merges a bunch of kvp files by a sort order
class merge_reader :  public kv_source
{
public:
	template<class It>
	merge_reader(
			const std::string& sort,
			const It& begin,
			const It& end,
			const std::string& begin_on,	
			const std::string& end_before,
			bool clean_break,
			const std::string& encoding);
	~merge_reader();
	bool read(std::string& key, std::string& value) override;

private:
	std::string m_begin_on;
	std::string m_end_before;
	bool m_clean_break;
	std::unique_ptr<sorter> m_sorter;
	std::vector<file_info_reader*> m_readers;  
	kv_merge* m_merged;
};

template<class It>
merge_reader::merge_reader(
		const std::string& sort,
		const It& begin,
		const It& end,
		const std::string& begin_on,
		const std::string& end_before,
		bool clean_break,
		const std::string& encoding)
	: m_begin_on(begin_on)
	, m_end_before(end_before)
	, m_clean_break(clean_break)
	, m_sorter(sorter_registry::get(sort, ""))
{
	// Build the sorter
	if (m_sorter == nullptr)
		throw io_exception("Unknown sorter");

	m_merged = new kv_merge(*m_sorter);     
	for(It cur = begin; cur != end; ++cur)
	{
		file_info_reader* reader = new file_info_reader(*cur, encoding);
		m_readers.push_back(reader);
		m_merged->add(reader);  
	}
}

#endif

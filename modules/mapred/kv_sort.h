
#ifndef __kv_sort_h__
#define __kv_sort_h__

#include "modules/io/keyvalue.h"
#include "modules/mapred/sorter.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/splitter.h"
#include <vector>
#include <queue>
#include <set>

class kv_sort : public kv_sink, public reset_kv_source
{
	struct kv_pair
	{
        	kv_pair(const std::string& key, const std::string& value);
        	std::string key;
        	std::string value;
	};
	typedef const kv_pair *pair_ptr;
	class ptr_compare
	{
	public: 
        	ptr_compare(std::shared_ptr<sorter> order);
        	bool operator()(pair_ptr a, pair_ptr b);
	private:
        	std::shared_ptr<sorter> m_order;
	};

	typedef std::vector<pair_ptr> collection_t;

public:
	struct param
	{
		param() {}
		param(const std::string& the_sorter, const std::string the_splitter)
			: sorter(the_sorter), splitter(the_splitter) {}
		
		TRANSFER_OBJECT
		{ 
			VERSION(0);
			FIELD(sorter, TF_ALLOW_NULL);
			FIELD(splitter, TF_ALLOW_NULL);
			FIELD(first_key, TF_ALLOW_NULL);
		}
		
		std::string sorter;
		std::string splitter;
		std::string first_key;
	};

	kv_sort(const std::string& serialized_params);
	virtual ~kv_sort();
	void write(const std::string& key, const std::string& value) override;
	void prep_read(); //potentially modify the underlying data before read ops start
	bool read(std::string& key, std::string& value) override;
	void reset() override;
	void clear();
	bool oversized(size_t goal_size) const { return (get_size() + 64 * get_num_records()) > goal_size; }
	bool legal_split(const std::string& key) const { return true; }
	size_t get_num_records() const { return m_pairs.size(); }
	size_t get_size() const { return m_records_size; }
	void set_file_info(file_info& fi) const;
	bool split_now(const std::string& key) const { return m_splitter ? (*m_splitter)(key) : false; }
	void update_split(const std::string& key) { if (m_splitter) m_splitter->set_initial_key(key); }

private:
	bool m_sorted;
	size_t m_records_size;
	std::shared_ptr<sorter> m_sorter;
	std::unique_ptr<splitter> m_splitter;
	ptr_compare m_order;
	collection_t m_pairs;
	collection_t::const_iterator m_current;
};

#endif

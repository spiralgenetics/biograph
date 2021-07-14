
#ifndef __kv_merge_h__
#define __kv_merge_h__

#include "modules/io/keyvalue.h"
#include "modules/mapred/sorter.h"
#include <vector>
#include <queue>

class file_info_reader;

class kv_merge : public kv_source
{
	struct kv_pair
	{
        	kv_pair(size_t part, const std::string& key, const std::string& value);
		size_t part;
        	std::string key;
        	std::string value;
	};
	typedef kv_pair *pair_ptr;
	class ptr_compare
	{
	public: 
        	ptr_compare(sorter& order);
        	bool operator()(pair_ptr a, pair_ptr b);
	private:
        	sorter& m_order;
	};

	typedef std::vector<pair_ptr> collection_t;
	typedef std::priority_queue<pair_ptr, collection_t, ptr_compare> queue_t;

public:
	kv_merge(sorter& order);
	~kv_merge();
	void add(file_info_reader* source);
	bool read(std::string& key, std::string& value) override;

private:
	queue_t m_queue;
	std::vector<file_info_reader*> m_sources;
};

#endif

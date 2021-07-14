
#ifndef __kv_summarize_h__
#define __kv_summarize_h__

#include "modules/io/transfer_object.h"
#include "modules/io/keyvalue.h"
#include "modules/mapred/sorter.h"
#include "modules/mapred/reducer.h"
#include "modules/mapred/manifest.h"
#include <vector>
#include <queue>
#include <set>

class kv_summarize : public kv_sink, public reset_kv_source
{
	class key_compare
	{
	public: 
        	key_compare(const sorter& order);
        	bool operator()(const std::string& a, const std::string& b);
	private:
        	const sorter& m_order;
	};

	typedef std::map<std::string, std::string, key_compare> collection_t;

public:
	class param
	{
	public:
		param() {}
		param(const std::string& in) { json_deserialize(*this, in); }

		TRANSFER_OBJECT
		{
			VERSION(0);
			FIELD(sort, TF_STRICT); 
			FIELD(reduce, TF_STRICT); 
			FIELD(reduce_param, TF_STRICT); 
		}
		std::string sort;
		std::string reduce;
		std::string reduce_param;
	};
	kv_summarize(const std::string& details);
	~kv_summarize();
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
	bool split_now(const std::string& key) const { return false; }
	void update_split(const std::string& key) {}

private:
	param m_params;
	size_t m_records_size;
	std::unique_ptr<sorter> m_sorter;
	std::unique_ptr<reducer> m_reducer;
	key_compare m_order;
	collection_t m_pairs;
	collection_t::const_iterator m_current;
};

#endif

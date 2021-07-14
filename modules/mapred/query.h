#ifndef __QUERY_H__
#define __QUERY_H__

#include "modules/io/keyvalue.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/mapred/manifest.h"

// typical usage:
// query qr;
// qr.find(sorted_data, myfirstKey, myLastKey);
// while( qr.read(k,v)) { //do something with (k,v) }

// query re-uses the same sorter that was used to order the manifest data.

class query : public kv_source
{
public:
	query();
	~query();

	bool read(std::string& key, std::string& value) override;

	// find the all key-value pairs included in [firstKey lastKey].
	// Each time find if called, the previous query results are deleted.
	// sorted_data must be globally sorted and firstKey < lastKey.
	void find(const manifest& sorted_data, const std::string& firstKey, const std::string& lastKey);

	template<typename TKeyType>
	void find_msgpack(const manifest& sorted_data, const TKeyType& firstKey, const TKeyType& lastKey)
	{
		find(sorted_data, msgpack_serialize(firstKey), msgpack_serialize(lastKey));
	}

	void clear(); //delete previous results

private:
	multi_reader<manifest::const_iterator>*	m_reader;
	std::string			m_firstKey,m_lastKey;
	bool				m_skippedToFirstKey;
	manifest			m_sorted_data;
};

#endif //__QUERY_H__

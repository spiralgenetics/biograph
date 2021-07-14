#include "modules/mapred/query.h"

query::query()
{
	m_reader = 0;
	m_skippedToFirstKey = false;
}

query::~query()
{
	clear();
}

//thanks to the pre-conditions verifications in query::find,
// we are guaranteed to have at least one value to read from a valid reader.
bool query::read(std::string& key, std::string& value)
{
	bool ret = false;
	if( m_reader )
	{
		std::string k,v;
		std::unique_ptr<sorter> s = sorter_registry::get(m_sorted_data.get_sort(), "");
		if( ! m_skippedToFirstKey )
		{// retrieve the first valid key-value pair
			do
			{//skip all until we get passed firstKey.
				ret = m_reader->read(k,v);
			}
			while( ret && s->lt( k, m_firstKey ) );

			if(s->gt(k, m_lastKey))
			{//this is the case when (firstKey,lastKey) both are not present in the dataset.
				clear();
				ret = false;
			}
			else
			{
				m_skippedToFirstKey = true;
				key = k; value = v;
				ret = true;
			}
		}
		else
		{// retrieve the rest of the key-value pairs
			ret = m_reader->read(k,v);

			if( !ret || s->gt( k, m_lastKey ) )
			{// we either reached the end or we passed the lastKey.
				std::string ik, iv;
				while(ret)  // Ignore it all
					ret = m_reader->read(ik,iv);
				clear();//we're through reading valid data.
			}
			else
			{
				key = k; value = v;
			}
		}
	}
	return ret;
}

void query::clear()
{
	if ( m_reader )
	{
		delete m_reader;
		m_reader = 0;
		m_firstKey = "";
		m_lastKey = "";
		m_skippedToFirstKey = false;
	}
}

void query::find(const manifest& sorted_data, const std::string& firstKey, const std::string& lastKey)
{
	clear(); //delete the previous query.

	manifest::const_iterator	it = sorted_data.begin(),
								end = sorted_data.end();
	if( it == end)
		return; //empty results

	//retrieve the sorter
	if( sorted_data.get_sort().empty() )
		return; // unsorted data

	std::unique_ptr<sorter> s = sorter_registry::get(sorted_data.get_sort(), "");

	if( s->gt(firstKey, lastKey) )
		return; //empty results;

	if( s->gt(it->first_key, lastKey) )
		return; // empty results;

	//find the first file that contains m_firstKey
	for ( ; it != end && s->lt(it->last_key, firstKey); it++);
	if(it == end)
		return; // empty results;

	manifest::const_iterator	first_file_info = it,
								end_file_info;

	//find the last file that may contains m_lastKey
	for ( ; it != end && !s->gt(it->last_key, lastKey); it++);

	// we want the one AFTER the last file containing lastKey, which
	// will be the 'end' point of the multi_reader.
	end_file_info = it == end ? end : ++it; //could be end, which is OK.

	std::string encoding = sorted_data.get_encoding();
	m_reader = make_multi_reader(first_file_info, end_file_info, encoding);

	//we'll need to remember the range when we read data from the reader.
	m_firstKey = firstKey;
	m_lastKey = lastKey;
	m_sorted_data = sorted_data;
}


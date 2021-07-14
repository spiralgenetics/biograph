#include "base/base.h"
#include "modules/bio_mapred/pair_reducer.h"

REGISTER_1(reducer, pair, const std::string&);


void pair_reducer::typed_start(read_id key)
{
	m_current_key = key;
	m_current_value.clear();
}


void pair_reducer::typed_add_value(const read_id& key, const unaligned_reads& value)
{
	if (! compare_read_id_stems(key, m_current_key)) {
		throw io_exception(
			boost::format("pair_reducer::add_value> key mismatch: current key is \"%1%\", new key is \"%2%\"")
				% m_current_key.pair_name % key.pair_name
		);
	}
	
	CHECK(value.size() == 1);
	m_current_value.push_back(value.at(0));
	CHECK(m_current_value.size() <= 2);
}


void pair_reducer::typed_end()
{
	output(m_current_key, m_current_value);
}


bool pair_reducer::compare_read_id_stems(const read_id& key1, const read_id& key2)
{
	std::string id1{key1.pair_name};
	std::string id2{key2.pair_name};
	
	if (id1.empty() || id2.empty()) {
		return id1 == id2;
	}
	
	if (*(id1.crbegin()) != '1' && *(id1.crbegin()) != '2') {
		return id1 == id2;
	}
	
	if (*(id2.crbegin()) != '1' && *(id2.crbegin()) != '2') {
		return id1 == id2;
	}
	
	id1.resize(id1.size() - 1);
	id2.resize(id2.size() - 1);

	return id1 == id2;
}

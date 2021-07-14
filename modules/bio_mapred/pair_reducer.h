#pragma once

#include "modules/mapred/reducer.h"
#include "modules/bio_base/unaligned_read.h"

class pair_reducer
	: public typed_reducer<pair_reducer, read_id, unaligned_reads,
							read_id, unaligned_reads>
{
public:
	pair_reducer(const std::string& params) {}

	virtual void typed_start(read_id key);
	virtual void typed_add_value(const read_id& key, const unaligned_reads& value);
	virtual void typed_end();

private:
	read_id m_current_key;
	unaligned_reads m_current_value;
	
	bool compare_read_id_stems(const read_id& key1, const read_id& key2);
};

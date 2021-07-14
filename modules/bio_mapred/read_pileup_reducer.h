
#pragma once

#include <memory>
#include "modules/bio_base/struct_var.h"
#include "modules/bio_base/var_info.h"
#include "modules/bio_base/pileup.h"
#include "modules/mapred/reducer.h"
#include "modules/mapred/query.h"
#include "modules/bio_base/reference.h"

struct read_pileup_params
{
	TRANSFER_OBJECT
        {
		VERSION(0);
		FIELD(reference);
		FIELD(min_depth);
		FIELD(var_infos);
        }

	std::string reference;
	uint32_t min_depth;
	manifest var_infos;
	void validate();
};

class read_pileup_reducer : public typed_reducer<read_pileup_reducer, struct_var_key, read_support, seq_position, struct_var>
{
public:
	read_pileup_reducer(const std::string& params);

	void typed_start(const struct_var_key& key);
	void typed_add_value(const struct_var_key& key, const read_support& value);
	void typed_end();

private:
	read_pileup_params m_params;
	uint32_t m_var_id;
	bool m_query_started;
	dna_sequence m_sequence;
	var_info m_var_info;
	query m_query;
	std::unique_ptr<reference> m_ref;
	std::unique_ptr<pileup> m_pileup;
};


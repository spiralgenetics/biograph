
#pragma once

#include "modules/io/transfer_object.h"

struct var_info
{
	TRANSFER_OBJECT { 
		VERSION(0);
		FIELD(s_ref);
		FIELD(e_ref);
		FIELD(avg_overlap);
		FIELD(s_flip);
		FIELD(e_flip);
		FIELD(is_ambig);
		FIELD(min_overlap);
	}
	
	uint32_t s_ref;
	uint32_t e_ref;
	float avg_overlap;
	bool s_flip;
	bool e_flip;
	bool is_ambig;
	uint8_t min_overlap;
};



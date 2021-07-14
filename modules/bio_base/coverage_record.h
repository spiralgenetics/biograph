
#pragma once

#include "modules/io/transfer_object.h"

struct coverage_record
{
	std::string read_name;
	size_t match_count;
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(read_name);
		FIELD(match_count);
	}
};

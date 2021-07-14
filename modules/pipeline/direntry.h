#pragma once

#include "modules/pipeline/datatype.h"
#include "modules/io/transfer_object.h"

class direntry
{
public:
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD_SPECIAL(COUCHDB_RESERVED, _id);
		FIELD_SPECIAL(COUCHDB_RESERVED, _rev);
		FIELD(url, TF_STRICT);
		FIELD(parent, TF_STRICT);
		FIELD(name, TF_STRICT);
		FIELD(user, TF_STRICT);
		FIELD(created, TF_STRICT);
		FIELD(directory, TF_STRICT);
		FIELD(type, TF_ALLOW_NULL);
		FIELD(size);
		FIELD(records);
		FIELD(reference);
		FIELD(sort_keys);
		FIELD(ref_name);
		FIELD(in_progress);
		FIELD(token);
	}

	std::string _id;
	std::string _rev;
	std::string url;
	std::string parent;
	std::string name;
	std::string user; // Empty if reference
	time_t created;
	bool directory;
	datatype_ref type;
	size_t size;
	size_t records;
	bool reference;
	std::vector<std::string> sort_keys;
	std::string ref_name;
	bool in_progress;
	std::string token;
};

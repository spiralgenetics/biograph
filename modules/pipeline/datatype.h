#pragma once

#include "modules/pipeline/restful_registry.h"
#include "modules/io/transfer_object.h"

struct datatype
{
	datatype(const std::string& id, const std::string& name, const std::string& description)
		: id(id)
		, name(name)
		, description(description)
	{}

	TRANSFER_OBJECT 
	{ 
		VERSION(0);
		FIELD(id, TF_STRICT); 
		FIELD(url, TF_STRICT); 
		FIELD(name, TF_STRICT); 
		FIELD(description, TF_STRICT); 
	}

	std::string id;
	std::string url;
	std::string name;
	std::string description;
};

typedef restful_registry<datatype> datatype_registry;
typedef datatype_registry::ref_type datatype_ref;

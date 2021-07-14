#pragma once

#include "modules/io/transfer_object.h"

// This is a variant structure that contains any additional data need to construct
// an importer or exporter.  The user sets the appropriate fields and the receiver
// is assumed to know how to unpack the data correctly.
struct ex_im_porter_data
{
	TRANSFER_OBJECT
	{
		VERSION(0);

		// The following optional fields are needed for importers/exporters.
		FIELD(ref_name);
		FIELD(skip_error);
		FIELD(start_key);
		FIELD(end_key);
	}

	std::string ref_name;
	bool skip_error;
	std::string start_key;
	std::string end_key;
};

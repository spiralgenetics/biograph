
#ifndef __seq_position_h__
#define __seq_position_h__

#include <boost/operators.hpp>

#include <string>
#include "modules/io/keyvalue.h"
#include "modules/io/transfer_object.h"

class seq_position
	: boost::less_than_comparable<seq_position>
	, boost::equality_comparable<seq_position>
{
public:
	seq_position();
	seq_position(int scaffold_id, unsigned long position);
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(scaffold_id);
		FIELD(position);
	};
	int scaffold_id;
	unsigned long position;

	bool valid() const { return scaffold_id != -1 && position != 0; }
	bool operator<(const seq_position& rhs) const;
	bool operator==(const seq_position& rhs) const;
	void bump_back(long dist);
};

#endif

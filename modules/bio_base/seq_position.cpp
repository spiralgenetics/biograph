
#include "modules/bio_base/seq_position.h"
#include <stdio.h>

seq_position::seq_position()
	: scaffold_id(-1)
	, position(0)
{}

seq_position::seq_position(int scaffold_id_, unsigned long position_)
	: scaffold_id(scaffold_id_)
	, position(position_)
{}

bool seq_position::operator<(const seq_position& rhs) const
{
	if (scaffold_id != rhs.scaffold_id)
		return scaffold_id < rhs.scaffold_id;
	return position < rhs.position;
}

bool seq_position::operator==(const seq_position& rhs) const
{
	return scaffold_id == rhs.scaffold_id && position == rhs.position;
}

void seq_position::bump_back(long dist)
{
  position = std::max(0L, long(position) - dist);
}


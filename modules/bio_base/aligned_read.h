#pragma once

#include "modules/bio_base/reference.h"
#include "modules/io/transfer_object.h"
#include "modules/bio_base/seq_position.h"

struct aligned_read
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(read_name);
		FIELD(flags);
		FIELD(ref_pos);
		FIELD(map_quality);
		FIELD(cigar);
		FIELD(mate_pos);
		FIELD(tlen);
		FIELD(seq);
		FIELD(qual);
		FIELD(read_group_id);
	};

	std::string read_name;
	int flags;
	seq_position ref_pos;
	int map_quality;
	std::string cigar;
	seq_position mate_pos;
	long tlen;
	std::string seq;
	std::string qual;
	std::string read_group_id;
};

bool parse_sam(const reference_assembly& ref, aligned_read& out, const std::string& sam_line);
std::string print_sam(const reference_assembly& ref, const aligned_read& in, bool use_supercontig_coords = false);

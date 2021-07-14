
#pragma once

#include "modules/io/transfer_object.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/unaligned_read.h"

struct corrected_read
{
	static const int no_match = 0;
	static const int start_match_fwd = 1;
	static const int start_match_rev = 2;
	static const int end_match_fwd = 3;
	static const int end_match_rev = 4;
	static const int invalid = 255;

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(quality);
		FIELD(sequence);
		FIELD(corrected);
		FIELD(match_type);
		FIELD(match_loc);
		FIELD(match_len);
		FIELD(trace_me);
		FIELD(aligned_pos);
	}

	std::string quality;          // Quality
	dna_sequence sequence;        // Original read sequence
	dna_sequence corrected;       // Corrected read sequence
	uint32_t match_loc = 0;       // Location of match (if any)
	uint8_t match_type = invalid; // Type of reference match (if any)
	uint8_t	match_len = 0;        // Length of match (if any)
	bool trace_me = false;
	seq_position aligned_pos;	  // sequence position according to aligner
};

typedef std::vector<corrected_read> corrected_reads;

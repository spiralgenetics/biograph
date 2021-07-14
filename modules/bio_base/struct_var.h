
#pragma once

#include "modules/io/transfer_object.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/seq_position.h"

struct struct_var_key
{
	struct_var_key() = default;
	struct_var_key(uint32_t the_variation_id, uint32_t the_read_id)
		: variation_id(the_variation_id)
		, read_id(the_read_id)
	{}

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(variation_id);
		FIELD(read_id);
	}

	bool operator<(const struct_var_key& rhs) const
	{
		if (variation_id < rhs.variation_id) return true;
		if (variation_id == rhs.variation_id) return read_id < rhs.read_id;
		return false;
	}
	bool operator==(const struct_var_key& rhs) const
	{
		return variation_id == rhs.variation_id &&
			read_id == rhs.read_id;
	}

	uint32_t variation_id = 0;
	uint32_t read_id = 0;
};

struct read_support
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(name);
		FIELD(original);
		FIELD(corrected);
		FIELD(quality);
		FIELD(pos);
		FIELD(flipped);
	}

	void flip();  // Flip

	std::string name;
	dna_sequence original;
	dna_sequence corrected;
	std::string quality;
	size_t pos;  // The position in the assembly
	bool flipped;  // Relative to original
};

struct struct_var
{
	enum ambside
	{
		amb_none,
		amb_left,
		amb_right,
	};

  // True if the assembled sequence doesn't match to something
  // nearby in the same extent and the same orientation.
	bool is_structural;
  // Points to the first base in reference that's outside the variant
  // region after trimming.  So, if the reference is "ABCDEFG" and the
  // variant is "ABCxyzEFG", ref_start would point to the "C" and
  // ref_end would point to the "E".
	seq_position ref_start;
	bool rev_start;
	seq_position ref_end;
	bool rev_end;

  // Per jeremy, probably not used.
	dna_sequence ref_seq;

  // Entire assembly, including some reference on both sides.
	dna_sequence assembled;

  // Within "assembled", the varying bases are within the half-open
  // range [var_start, var_end).
	size_t var_start;  // Within assembled sequence
	size_t var_end;    // Within assembled sequence

  // Minimum coverage depth of variant, including a base in reference
  // on each side.
	size_t depth;
  // Unique identifier for this assembly, passed through to generated files.
	uint32_t var_id;

	bool flipped;  // Relative to supporting reads
	bool is_ambig;

	double avg_depth;
	uint8_t min_overlap;
	double avg_overlap;
	std::vector<size_t> assembly_depth;
	std::vector<size_t> assembly_fwd;
	std::vector<size_t> assembly_tot_qual;

  // If there are ares of 0 coverage, has_holes is true.  An variant
  // should not be called in this case.
	bool has_holes; // Did the assembly have holes?

  // If the cost gets too hi during a*, this flag is set and a larger
  // variant is returned.
	bool align_failed; // Did we fail to align (ie, possible cluster)

	int sub_id; // Subassembly id
	std::string filter;
	size_t ambiguous_side = amb_none; // is one of ambside
	size_t ambiguous_count = 0;
	std::string transpose;  // If right side hits ALU db, scaffold_id of right hand ALU db hit
	double simple_alignment_score;

	void flip();
	void canonicalize(); // Flip as needed so start < end

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(is_structural);
		FIELD(ref_start);
		FIELD(rev_start);
		FIELD(ref_end);
		FIELD(rev_end);
		FIELD(ref_seq);
		FIELD(assembled);
		FIELD(var_start);
		FIELD(var_end);
		FIELD(depth);
		FIELD(var_id);
		FIELD(flipped);
		FIELD(is_ambig);
		FIELD(avg_depth);
		FIELD(min_overlap);
		FIELD(avg_overlap);
		FIELD(assembly_depth);
		FIELD(has_holes);
		FIELD(align_failed);
		FIELD(sub_id);
		FIELD(filter);
		FIELD(ambiguous_side);
		FIELD(ambiguous_count);
		FIELD(transpose);
		FIELD(simple_alignment_score);
		FIELD(assembly_fwd);
		FIELD(assembly_tot_qual);
	}
};

typedef std::vector<struct_var> struct_vars;

// return true if pos is 100bp left of the beginning of the scaffold
// and 100bp right of its end
bool safe_range(size_t pos, size_t tot_size);


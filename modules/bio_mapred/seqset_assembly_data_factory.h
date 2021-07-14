#pragma once

#include "modules/io/version.h"

const product_version k_assembly_data_version{"1.0.0"};


// The plain old data portion of the assembly data, i.e. everything except the
// bitcount buffers and reverse complement vector.
struct seqset_assembly_pod
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(version, TF_STRICT);
		FIELD(min_overlap, TF_STRICT);
		FIELD(max_overlap_count, TF_STRICT);
		FIELD(read_length, TF_STRICT);
	}

	product_version version;
	int min_overlap; // Min number of bases that must overlap to proceed with assembly.
	int max_overlap_count; // Max number of overlaps before we give up looking for more.
	unsigned read_length;

	seqset_assembly_pod() : version(k_assembly_data_version) {}
	seqset_assembly_pod(
		int the_min_overlap
		, int the_max_overlap_count
		, unsigned the_read_length
	) : version{k_assembly_data_version}
	, min_overlap(the_min_overlap)
	, max_overlap_count(the_max_overlap_count)
	, read_length(the_read_length)
	{}
};

inline bool operator==(const seqset_assembly_pod& lhs, const seqset_assembly_pod& rhs)
{
	return lhs.version == rhs.version
		&& lhs.min_overlap == rhs.min_overlap
		&& lhs.max_overlap_count == rhs.max_overlap_count
		&& lhs.read_length == rhs.read_length
	;
}

class seqset;
class seqset_bitmap_base;

// Contains additional data structures needed to efficiently assemble reads from a seqset.
// This struct will be built in memory and then copied into the memmap because we don't
// know in advance how big the bitcount buffers will have to be until they are actually
// built.
struct seqset_assembly_data_factory
{
public:
	seqset_assembly_data_factory(const seqset& seqset, int min_overlap, int max_overlap_count
						, const seqset_bitmap_base& the_bitmap)
		: m_seqset(seqset)
		, m_bitmap(the_bitmap)
		, m_pod{min_overlap, max_overlap_count, m_seqset.read_len()}
	{ build_data_structures(); }

	const seqset& m_seqset;
	const seqset_bitmap_base& m_bitmap;
	seqset_assembly_pod m_pod;

	std::vector<char> m_read_bitcount_buffer;
	// A bit for every seqset entry, set if the entry is a full-length read.
	std::unique_ptr<bitcount> m_read_bitcount;

	// Table of reverse complements for every read.  The index (read ID) is given by find_count in the
	// read bitcount, i.e. it's the number of the read in the seqset listing them in entry order.
	// The value is the read ID of the read reverse complement
	std::vector<uint32_t> m_read_rcs;

	std::vector<char> m_nonunique_overlap_buffer;
	// A bitcount for every read with a bit reset (0) when there is a unique overlap
	// in both the right and left direction (white node).  A 1 bit (black node)
	// can mean no overlaps or multiple overlaps in one or both directions.
	std::unique_ptr<bitcount> m_nonunique_overlap;

	// A table of the black nodes that can be reached via overlaps from a given black
	// node.  The data table contains all of the black node targets concatenated and
	// the offset table is the index into the data table for a give node.   The offset
	// is terminated by the size of the data table as an extra, final entry.  The black
	// nodes are stored by node_id, their count in the nonunique bitmap.
	std::vector<uint32_t> m_left_black_nodes_data;
	std::vector<uint32_t> m_left_black_nodes_offsets;

	// An accessor function for m_left_black_nodes_data.  Pass in a node ID, get a
	// vector of terminal black nodes for all the node's overlaps.
	std::vector<uint32_t> get_reachable_black_nodes(uint32_t node_id);

	// Take a hash table returned from seqset_range::find_overlap_reads and follow the
	// overlaps in any white nodes to black nodes.  Returns a vector of read IDs (not
	// entries!) for the black nodes.
	std::vector<uint32_t> follow_to_black(overlaps_t& results) const;

	// Follow the overlaps from a single node from a seqset entry until a black node is
	// reached.  If the passed-in entry is black, the function returns its input.
	// Returns the read ID of the final black node.
	uint32_t follow_one_to_black(uint64_t entry) const;

	bool is_node_white(uint64_t node_entry) const;

private:
	void build_data_structures();
	void build_read_bitmap();
	void build_rc_table();
	void build_unique_overlaps();
	void build_unique_assemblies();
};

#pragma once

#include <string>

#include "modules/io/bitcount.h"
#include "modules/bio_mapred/seqset_assembly_data_factory.h"

struct alignas(8) seqset_assembly_offsets
{
	size_t m_read_bitcount_offset;
	size_t m_unique_overlap_offset;
	size_t m_rc_array_offset;
	size_t m_left_black_offsets_offset;
	size_t m_left_black_data_offset;
	size_t m_serialized_pod_offset;

	seqset_assembly_offsets(const void * seqset_assembly_offsets_ptr)
	{
		typed_memcpy(this, static_cast<const seqset_assembly_offsets*>(seqset_assembly_offsets_ptr));
	}

	seqset_assembly_offsets(
		size_t entry_count
		, size_t read_count
		, size_t black_node_count
		, size_t black_node_data_size
	) : m_read_bitcount_offset(sizeof(*this))
	, m_unique_overlap_offset(m_read_bitcount_offset + bitcount::compute_size(entry_count))
	, m_rc_array_offset(m_unique_overlap_offset + bitcount::compute_size(read_count))
	, m_left_black_offsets_offset(m_rc_array_offset + read_count * sizeof(uint32_t))
	, m_left_black_data_offset(m_left_black_offsets_offset + black_node_count * sizeof(uint32_t))
	, m_serialized_pod_offset(m_left_black_data_offset + black_node_data_size * sizeof(uint32_t))
	{}

	size_t get_last_offset() const { return m_serialized_pod_offset; }
};

class seqset;
class seqset_bitmap_base;
struct seqset_assembly_data
{
	seqset_assembly_data(
		const std::string& mmap_file_path
		, const seqset& seqset
		, const seqset_bitmap_base& the_bitmap
	);

	// This function builds the data structures needed to assemble from a seqset, and then builds
	// a memmap to hold the data and writes it all to disk. The written mmap file can then be
	// used to construct a seqset_assembly_data object.
	static void build_seqset_assembly_data(
		const std::string& file_path
		, const seqset& seqset
		, const seqset_bitmap_base& the_bitmap
		, int min_overlap
		, int max_overlap_count
	);

	// Same as follow_one_to_black except that this function returns the assembled sequence
	// and a vector of overlaps.
	using assembly_info_t = std::pair<dna_sequence, std::vector<uint8_t>>;
	assembly_info_t assemble_to_black(uint64_t entry) const;

	bool is_node_white(uint64_t node_entry) const;
	uint32_t get_read_id(uint64_t node_entry) const { return m_read_bitcount.count(node_entry); }
	uint32_t get_black_node_id(uint64_t node_entry) const { return m_nonunique_overlap.count(get_read_id(node_entry)); }
	uint32_t read_to_black_node_id(uint32_t read_id) const { return m_nonunique_overlap.count(read_id); }
	std::vector<uint32_t> get_assembled_node_ids(uint32_t node_id) const
	{
		return std::vector<uint32_t>(m_left_black_nodes_offset_array + node_id, m_left_black_nodes_offset_array + node_id + 1);
	}

	const mmap_buffer m_mmap;
	const seqset& m_seqset;
	const seqset_bitmap_base& m_bitmap;
	const seqset_assembly_offsets m_offsets;
	const bitcount m_read_bitcount;
	const bitcount m_nonunique_overlap;
	const uint32_t* m_rc_array; // Has size m_read_bitcount.total_bits()
	const uint32_t* m_left_black_nodes_offset_array; // Has size m_nonunique_overlap.total_bits() + 1
	const uint32_t* m_left_black_nodes_data_array; // Size is last element of m_left_black_nodes_offset_array
	const seqset_assembly_pod m_pod;
};

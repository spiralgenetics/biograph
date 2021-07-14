#include <vector>
#include <memory>
#include <cstring>

#include "modules/io/mmap_buffer.h"
#include "modules/io/bitcount.h"
#include "modules/io/utils.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_mapred/seqset_assembly_data_factory.h"
#include "modules/bio_mapred/seqset_assembly_data.h"

seqset_assembly_data::seqset_assembly_data(
	const std::string& mmap_file_path
	, const seqset& seqset
	, const seqset_bitmap_base& the_bitmap
) : m_mmap{mmap_file_path}, m_seqset(seqset), m_bitmap(the_bitmap)
	, m_offsets(m_mmap.buffer())
	, m_read_bitcount{m_mmap.buffer() + m_offsets.m_read_bitcount_offset
		, seqset.size()}
	, m_nonunique_overlap{m_mmap.buffer() + m_offsets.m_unique_overlap_offset
		, m_read_bitcount.total_bits()}
	, m_rc_array(reinterpret_cast<const uint32_t*>(m_mmap.buffer() + m_offsets.m_rc_array_offset))
	, m_left_black_nodes_offset_array(reinterpret_cast<const uint32_t*>(m_mmap.buffer()
		+ m_offsets.m_left_black_offsets_offset))
	, m_left_black_nodes_data_array(reinterpret_cast<const uint32_t*>(m_mmap.buffer()
		+ m_offsets.m_left_black_data_offset))
{
	std::string serialized_pod{m_mmap.buffer() + m_offsets.m_serialized_pod_offset
		, m_mmap.size() - m_offsets.m_serialized_pod_offset};
	msgpack_deserialize(m_pod, serialized_pod);
}

void seqset_assembly_data::build_seqset_assembly_data(
	const std::string& file_path
	, const seqset& seqset
	, const seqset_bitmap_base& the_bitmap
	, int min_overlap
	, int max_overlap_count
)
{
	seqset_assembly_data_factory memory_assembly_data{
		seqset, min_overlap, max_overlap_count, the_bitmap};

	seqset_assembly_offsets memmap_offsets(
		memory_assembly_data.m_read_bitcount->size() // Number of seqset entries
		, memory_assembly_data.m_nonunique_overlap->size() // Number of reads
		, memory_assembly_data.m_nonunique_overlap->total_bits() + 1 // Number of black nodes
		, memory_assembly_data.m_left_black_nodes_data.size()
	);

	std::string serialized_pod{ msgpack_serialize(memory_assembly_data.m_pod) };
	size_t assembly_data_size = memmap_offsets.get_last_offset() + serialized_pod.size();

	mmap_buffer assembly_data_mmap{file_path, assembly_data_size};

	typed_memcpy(assembly_data_mmap.buffer(), &memmap_offsets);
	vector_memcpy(assembly_data_mmap.buffer() + memmap_offsets.m_read_bitcount_offset
		, memory_assembly_data.m_read_bitcount_buffer);
	vector_memcpy(assembly_data_mmap.buffer() + memmap_offsets.m_unique_overlap_offset
		, memory_assembly_data.m_nonunique_overlap_buffer);
	vector_memcpy(assembly_data_mmap.buffer() + memmap_offsets.m_rc_array_offset
		, memory_assembly_data.m_read_rcs);
	vector_memcpy(assembly_data_mmap.buffer() + memmap_offsets.m_left_black_offsets_offset
		, memory_assembly_data.m_left_black_nodes_offsets);
	vector_memcpy(assembly_data_mmap.buffer() + memmap_offsets.m_left_black_data_offset
		, memory_assembly_data.m_left_black_nodes_data);
	vector_memcpy(assembly_data_mmap.buffer() + memmap_offsets.m_serialized_pod_offset
		, serialized_pod);
}

seqset_assembly_data::assembly_info_t seqset_assembly_data::assemble_to_black(uint64_t node_entry) const
{
	seqset_assembly_data::assembly_info_t assembly_info;

	seqset_range node_context = m_seqset.ctx_entry(node_entry);
	std::get<0>(assembly_info) = node_context.sequence();
	overlaps_t node_overlaps;
	std::vector<uint8_t> overlaps;
	while (is_node_white(node_entry)) {
		node_overlaps.clear();
		node_context.find_overlap_reads(node_overlaps, m_pod.max_overlap_count, m_pod.min_overlap, m_bitmap);
		overlaps.push_back(node_overlaps.cbegin()->second);
		node_entry = node_overlaps.cbegin()->first;
		node_context = m_seqset.ctx_entry(node_entry);
		std::get<0>(assembly_info) = node_context.sequence()
			.subseq(0, m_pod.read_length - overlaps.back()) + std::get<0>(assembly_info);
	}

	return assembly_info;
}

bool seqset_assembly_data::is_node_white(uint64_t node_entry) const
{
	uint32_t node_read_id = m_read_bitcount.count(node_entry);
	return ! m_nonunique_overlap.get(node_read_id);
}

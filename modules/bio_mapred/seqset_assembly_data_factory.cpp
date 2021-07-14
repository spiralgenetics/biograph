#include "base/base.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_mapred/seqset_assembly_data_factory.h"

#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_thread_num() { return 0; }
#endif

void seqset_assembly_data_factory::build_data_structures()
{
	SPLOG("Building assembly data with read length = %u, min_overlap = %d and max_overlap_count = %d"
		, m_pod.read_length, m_pod.min_overlap, m_pod.max_overlap_count);
	build_read_bitmap();
	build_rc_table();
	build_unique_overlaps();
	build_unique_assemblies();
}

void seqset_assembly_data_factory::build_read_bitmap()
{
	SPLOG("Populating pop_front cache...");
	m_seqset.populate_pop_front_cache();

	m_read_bitcount_buffer.resize(bitcount::compute_size(m_seqset.size()));
	m_read_bitcount = make_unique<bitcount>(m_read_bitcount_buffer.data(), m_seqset.size());
	m_read_bitcount->init();

	SPLOG("Generating sequence table from %lu seqset entries with min_overlap = %d",
		m_seqset.size(), m_pod.min_overlap);
	for (auto i = 0UL; i < m_seqset.size(); i++) {
		seqset_range the_context{ m_seqset.ctx_entry(i) };
		if (the_context.size() == m_pod.read_length && m_bitmap.get_bit(i)) {
			m_read_bitcount->set(i, true);
		} else {
			m_read_bitcount->set(i, false);
		}
	}
	m_read_bitcount->finalize();
	SPLOG("Read bitcount complete with %lu entries and %lu reads.", m_read_bitcount->size(), m_read_bitcount->total_bits());
}

void seqset_assembly_data_factory::build_rc_table()
{
	m_read_rcs.resize(m_read_bitcount->total_bits());
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for (auto i = 0UL; i < m_seqset.size(); i++) {
		if (! m_read_bitcount->get(i)) { continue; }
		seqset_range the_context{ m_seqset.ctx_entry(i) };
		seqset_range rc_context{ m_seqset.find(the_context.sequence().rev_comp()) };
		CHECK(rc_context.valid());
		CHECK(rc_context.end() - rc_context.begin() == 1);
		m_read_rcs[m_read_bitcount->count(the_context.begin())] = m_read_bitcount->count(rc_context.begin());
	}
	SPLOG("Built reverse complement table with %lu entries", m_read_rcs.size());
}

void seqset_assembly_data_factory::build_unique_overlaps()
{
	SPLOG("Building unique overlap readmap");
	m_nonunique_overlap_buffer.resize(bitcount::compute_size(m_read_bitcount->total_bits()));
	m_nonunique_overlap = make_unique<bitcount>(m_nonunique_overlap_buffer.data(), m_read_bitcount->total_bits());
	m_nonunique_overlap->init();

	SPLOG("Looking for left unique overlaps.");
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for (auto i = 0UL; i < m_seqset.size(); i++) {
		if (! m_read_bitcount->get(i)) { continue; }
		overlaps_t overlap_result;
		seqset_range overlap_context{ m_seqset.ctx_entry(i) };
		bool is_left_unique_overlap =
			(overlap_context.find_overlap_reads(overlap_result, 1, m_pod.min_overlap, m_bitmap))
			&& (! overlap_result.empty());
		m_nonunique_overlap->set(m_read_bitcount->count(i), ! is_left_unique_overlap);
	}
	m_nonunique_overlap->finalize();
	SPLOG("Leftwards non-unique overlap readmap has %lu entries and %lu non-unique left overlaps."
		, m_nonunique_overlap->size()
		, m_nonunique_overlap->total_bits());

	SPLOG("Looking for right unique overlaps.");
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for (auto i = 0UL; i < m_seqset.size(); i++) {
		if (! m_read_bitcount->get(i)) { continue; }
		auto read_id = m_read_bitcount->count(i);
		if (m_nonunique_overlap->get(read_id)) { continue; }
		if (m_nonunique_overlap->get(m_read_rcs[read_id])) {
			m_nonunique_overlap->set(read_id, true);
		}
	}

	m_nonunique_overlap->finalize();
	SPLOG("Non-unique overlap readmap complete with %lu entries and %lu non-unique overlaps in both directions."
		, m_nonunique_overlap->size()
		, m_nonunique_overlap->total_bits());
}

void seqset_assembly_data_factory::build_unique_assemblies()
{
	std::vector<std::vector<uint32_t>> left_black_node_assemblies(m_nonunique_overlap->total_bits());
	SPLOG("Resized black node offsets for %lu left black nodes.", m_left_black_nodes_offsets.size());
	int thread0_count = 0;
	std::atomic_size_t total_target_node_count{0};
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for (auto i = 0UL; i < m_nonunique_overlap->size(); i++) {
		if (! m_nonunique_overlap->get(i)) { continue; }
		if (omp_get_thread_num() == 0 && (++thread0_count % 100000 == 0)) {
			SPLOG("Processing node %lu on %d iterations", i, thread0_count);
		}
		auto read_id = m_read_bitcount->find_count(i);
		seqset_range read_context{ m_seqset.ctx_entry(read_id) };
		overlaps_t left_overlaps;
		read_context.find_overlap_reads(left_overlaps, m_pod.max_overlap_count, m_pod.min_overlap, m_bitmap);
		left_black_node_assemblies[m_nonunique_overlap->count(i)] = follow_to_black(left_overlaps);
		total_target_node_count.fetch_add(left_overlaps.size());
	}

	SPLOG("Flattening left black node assemblies.");
	m_left_black_nodes_offsets.reserve(m_nonunique_overlap->total_bits() + 1);
	m_left_black_nodes_offsets.push_back(0U);
	m_left_black_nodes_data.reserve(total_target_node_count);
	for (auto& assembled_nodes : left_black_node_assemblies) {
		m_left_black_nodes_offsets.push_back(m_left_black_nodes_offsets.back() + assembled_nodes.size());
		std::move(assembled_nodes.begin(), assembled_nodes.end(), std::back_inserter(m_left_black_nodes_data));
	}
	SPLOG("%lu black nodes assembled.", m_left_black_nodes_offsets.size() - 1);
}

std::vector<uint32_t> seqset_assembly_data_factory::get_reachable_black_nodes(uint32_t node_id)
{
	return std::vector<uint32_t>(
		m_left_black_nodes_data.begin() + m_left_black_nodes_offsets[node_id]
		, m_left_black_nodes_data.begin() + m_left_black_nodes_offsets[node_id + 1]
	);
}

std::vector<uint32_t> seqset_assembly_data_factory::follow_to_black(overlaps_t& results) const
{
	using overlap_kvp = std::pair<uint64_t, uint8_t>;

	std::vector<uint32_t> black_nodes;
	black_nodes.reserve(results.size());
	std::transform(
		results.begin()
		, results.end()
		, std::back_inserter(black_nodes)
		, [this](const overlap_kvp& the_overlap_kvp) { return follow_one_to_black(the_overlap_kvp.first); }
	);

	return black_nodes;
}

uint32_t seqset_assembly_data_factory::follow_one_to_black(uint64_t node_entry) const
{
	while (is_node_white(node_entry)) {
		seqset_range node_context = m_seqset.ctx_entry(node_entry);
		overlaps_t node_overlaps;
		node_context.find_overlap_reads(node_overlaps, m_pod.max_overlap_count, m_pod.min_overlap, m_bitmap);
		node_entry = node_overlaps.cbegin()->first;
	}

	return m_nonunique_overlap->count(m_read_bitcount->count(node_entry));
}

bool seqset_assembly_data_factory::is_node_white(uint64_t node_entry) const
{
	uint32_t node_read_id = m_read_bitcount->count(node_entry);
	return ! m_nonunique_overlap->get(node_read_id);
}

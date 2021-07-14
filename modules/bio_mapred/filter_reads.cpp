
#include "base/base.h"
#include "modules/bio_mapred/filter_reads.h"
#include "modules/io/json_transfer.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/kmer.h"

REGISTER_1(mapper, filter_reads, const std::string&);

filter_reads_mapper::filter_reads_mapper(const std::string& params)
	: m_tot_mapped(0)
	, m_tot_filtered(0)
{
	json_deserialize(m_params, params);
	m_params.validate();
}

void filter_reads_mapper::setup()
{
	m_kdb = std::unique_ptr<kmer_set>(new kmer_set(m_params.kmer_db, [&](double) {m_watchdog();}));
}

void filter_reads_mapper::typed_map(const std::string& key, const corrected_reads& value)
{
	corrected_reads outs;
	for(size_t i = 0; i < value.size(); i++)
	{
		corrected_read out = value[i];
		if (map_one_read(key, out))
			outs.push_back(out);
	}

	if (!outs.empty() )
	{
		output(key, outs);
	}
}

bool filter_reads_mapper::map_one_read(const std::string& key, corrected_read& r)
{
	if (r.corrected.size() > m_mismatch_count.size()) {
		m_mismatch_count.resize(r.corrected.size());
	}

	unsigned mismatches = 0;
	size_t kmer_size = m_kdb->kmer_size();
	for(size_t i = 0; i < r.corrected.size() - kmer_size; i++) {
		kmer_t k = make_kmer(r.corrected.begin() + i, kmer_size);
		k = canonicalize(k, kmer_size);
		if (!m_kdb->count(k)) {
			mismatches++;
		}
	}

	CHECK(mismatches < m_mismatch_count.size());
	m_mismatch_count[mismatches]++;
	m_tot_mapped++;
	if (mismatches) {
		r.trace_me = true;
		m_tot_filtered++;
	}
	return !m_params.hard_filter || (mismatches > 0);
}

void filter_reads_mapper::install_metadata(meta::data& metadata)
{
	SPLOG("Goodbye cruel world: m_tot_mapped = %zu, m_tot_filtered = %zu",
		m_tot_mapped, m_tot_filtered);

	// Truncate trailing zeroes from distribution vector.
	auto last_non_zero_iter = std::find_if(m_mismatch_count.rbegin(), m_mismatch_count.rend(),
		[](unsigned element) { return element != 0; }
	).base();
	if (last_non_zero_iter != m_mismatch_count.end() && last_non_zero_iter != m_mismatch_count.begin()) {
		m_mismatch_count.erase(++last_non_zero_iter, m_mismatch_count.end());
	}

	metadata.set(meta::ns::readonly, "tagged_reads_count", m_tot_filtered);
	metadata.set(meta::ns::readonly, "filtered_read_dist", m_mismatch_count);
}

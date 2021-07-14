
#include "modules/bio_mapred/filter_assemblies.h"
#include "modules/io/json_transfer.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/kmer.h"

REGISTER_1(mapper, filter_assemblies, const std::string&);

filter_assemblies_mapper::filter_assemblies_mapper(const std::string& params)
{
	json_deserialize(m_params, params);
	m_params.validate();
}

filter_assemblies_mapper::~filter_assemblies_mapper()
{}

void filter_assemblies_mapper::setup()
{
	SPLOG("filter_assemblies_mapper::setup> Loading kmers")
	m_kdb = make_unique<kmer_set>(m_params.kmer_db, [&](double) { m_watchdog(); });
	SPLOG("filter_assemblies_mapper::setup> Kmers loaded")
}

void filter_assemblies_mapper::typed_map(const seq_position& key, const struct_var& sv)
{
	size_t interesting = 0;
	size_t kmer_size = m_kdb->kmer_size();
	for(size_t i = sv.var_start - kmer_size + 1; i < sv.var_end; i++) {
		kmer_t k = make_kmer(sv.assembled.begin() + i, kmer_size);
		k = canonicalize(k, kmer_size);
		if (!m_kdb->count(k)) {
			interesting++;
		}
	}
	m_mapped_count++;

	if (interesting >= kmer_size - 5) {
		m_filtered_count++;
		output(key, sv);
	}
	
	if (interesting > 100) {
		interesting = 100;
	}
	if (interesting >= m_mismatch_distribution.size()) {
		m_mismatch_distribution.resize(interesting + 1);
	}
	m_mismatch_distribution[interesting]++;
}

void filter_assemblies_mapper::install_metadata(meta::data& metadata)
{
	SPLOG("filter_assemblies_mapper::install_metadata> Mapped %zu reads, passed %zu through filter.",
		m_mapped_count, m_filtered_count); 

	// Truncate trailing zeroes from distribution vector.
	auto last_non_zero_iter = std::find_if(m_mismatch_distribution.rbegin(), m_mismatch_distribution.rend(),
		[](unsigned element) { return element != 0; }
	).base();
	if (last_non_zero_iter != m_mismatch_distribution.end() && last_non_zero_iter != m_mismatch_distribution.begin()) {
		m_mismatch_distribution.erase(++last_non_zero_iter, m_mismatch_distribution.end());
	}

	metadata.set(meta::ns::readonly, "tagged_assembly_count", m_filtered_count);
	metadata.set(meta::ns::readonly, "filtered_assembly_dist", m_mismatch_distribution);
}

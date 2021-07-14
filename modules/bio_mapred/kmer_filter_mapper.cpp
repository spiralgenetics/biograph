
#include "modules/mapred/manifest.h"
#include "modules/bio_mapred/kmer_filter_mapper.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"

REGISTER_1(mapper, kmer_filter, const std::string&);

void kmer_filter_params::validate()
{
	SPLOG_P(LOG_DEBUG, "kmer_filter_params::validate> min_count: %lu", min_count);
}

kmer_filter_mapper::kmer_filter_mapper(const std::string& params)
	: m_params(inline_json_deserialize<kmer_filter_params>(params))
	, m_overrep(m_params.kmer_size)
{
	m_params.validate();
}

void kmer_filter_mapper::setup()
{
	manifest_reader mr(m_params.overrep);
	kmer_t k;
	kcount_pair v;
	SPLOG_P(LOG_DEBUG, "kmer_filter_mapper::setup> Reading overrep data");
	while(mr.read_msgpack(k, v)) {
		m_overrep.add_overrep(overrep_t(k, v.fwd + v.rev));
		//SPLOG("%s: %u,%u", dna_sequence(k, m_params.kmer_size).as_string().c_str(), v.fwd, v.rev);
		m_watchdog();
	}
	SPLOG_P(LOG_DEBUG, "kmer_filter_mapper::setup> Done reading overrep data, count = %lu", m_overrep.size());
}

void kmer_filter_mapper::typed_map(kmer_t kmer, const kcount_pair& count)
{
	overrep_t o;
	if (m_overrep.find_near(kmer, o)) {
		uint32_t min_c = std::min(count.fwd, count.rev);
		uint32_t max_c = std::max(count.fwd, count.rev);
		if (min_c < o.second * m_params.rnd_err_thresh &&
			max_c < o.second * m_params.sys_err_thresh) {
			//SPLOG("Dropping: kmer: %s (%d, %d), overrep: %s (%d)",
			//	kmer_str(kmer, ks).c_str(), count.fwd, count.rev,
			//	kmer_str(o.first, ks).c_str(), o.second);
			return;
		}
	}
	if (count.fwd + count.rev > m_params.min_count)
		output(kmer, count);
}



#include "modules/bio_base/reference.h"
#include "modules/bio_mapred/correct_reads_mapper.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/bio_base/fast_read_correct.h"
#include "modules/mapred/metadata.h"

typedef correct_reads_mapper correct_reads_dual_mapper;
REGISTER_1(dual_mapper, correct_reads, const std::string&);

meta::merge::init merge_read_count(meta::merge::register_fn("corrected_read_count",
	meta::merge::sum
));

meta::merge::init merge_read_bases(meta::merge::register_fn("corrected_read_bases",
	meta::merge::sum
));

meta::merge::init merge_failed_count(meta::merge::register_fn("failed_correction_count",
	meta::merge::sum
));

meta::merge::init merge_base_dist(meta::merge::register_fn("corrected_base_dist", [] (
		const meta::merge::params& params)
	{
		std::vector<uint64_t> value1;
		std::vector<uint64_t> value2;
		json_unwrap(value1, params.value1);
		json_unwrap(value2, params.value2);

		if (value1.size() < value2.size()) {
			value1.swap(value2);
		}

		std::vector<uint64_t> result;
		result.reserve(value1.size());

		std::transform(
			value1.begin(),
			value1.begin() + value2.size(),
			value2.begin(),
			std::back_inserter(result), std::plus<uint64_t>()
		);

		return json_wrap(result);
	}
));

correct_reads_mapper::correct_reads_mapper(const std::string& params)
{
	json_deserialize(m_params, params);
	m_params.validate();
	SPLOG_P(LOG_DEBUG, "correct_reads_mapper> trim: %d", (int) m_params.trim);
	SPLOG_P(LOG_DEBUG, "correct_reads_mapper> max_quality_cost=%.0f", m_params.max_quality_cost);
}

task_requirements correct_reads_mapper::get_requirements()
{
	// SPLOG_P(LOG_DEBUG, "correct_reads_mapper::get_requirements>");
	return task_requirements {
		.profile = "normal",
		.cpu_minutes = 3,
	};
}

void correct_reads_mapper::setup()
{
	// SPLOG_P(LOG_DEBUG, "correct_reads_mapper::setup>");
	m_kdb = make_unique<kmer_set>(m_params.kmer_db, [&](double progress) {
		m_watchdog();
	});
}

void correct_reads_mapper::typed_map(const read_id& key, const unaligned_reads& value)
{
	// SPLOG_P(LOG_DEBUG, "correct_reads_mapper::typed_map>");
	corrected_reads outs;
	for(size_t i = 0; i < value.size(); i++)
	{
		corrected_read out;
		if (map_one_read(out, key, value[i]))
			outs.push_back(out);
	}

	if (! outs.empty() )
	{
		output1(key.pair_name, outs);
	}
}

bool correct_reads_mapper::map_one_read(corrected_read& v, const read_id& id, const unaligned_read& r)
{
	// SPLOG_P(LOG_DEBUG, "correct_reads_mapper::map_one_read>");
	v.quality = r.quality;
	v.sequence = dna_sequence(r.sequence);
	v.aligned_pos = r.ref_loc;

	if (m_params.skip_snps && (r.mismatches >= 0 && r.mismatches <= 2)) {
		coverage_record cr;
		cr.read_name = id.pair_name;
		cr.match_count = 1;
		output2(r.ref_loc, cr);
		return false;
	}

	if (m_params.trim != 0) {
		if (r.sequence.size() < m_params.trim) {
			return false;
		}
		v.quality = v.quality.substr(0, m_params.trim);
		v.sequence = v.sequence.subseq(0, m_params.trim);
	}

	if (r.sequence.size() < m_kdb->kmer_size())
		return false;

   if (m_params.frc_max_corrections)  {
     frc_params params;
     params.max_corrections = m_params.frc_max_corrections;
     params.min_good_run = m_params.frc_min_good_run;
     params.kmer_size = m_kdb->kmer_size();
     params.kmer_lookup_f = [&](kmer_t kmer, frc_kmer* ki) -> bool {
       kmer_t canon = canonicalize(kmer, m_kdb->kmer_size(), ki->flipped);
       auto index = m_kdb->find_table_index(canon);
       if (index == kmer_set::k_not_present) {
         return false;
       }
       ki->index = index;
       return true;
     };
     frc_output res = fast_read_correct(r.sequence, params);
     unsigned needed_good_bases = m_params.trim_after_portion * v.sequence.size();
     if (res.corrected.size() < needed_good_bases) {
       m_stats.failed_correction_count++;
       return false;
     }

     v.corrected = std::move(res.corrected);
     if (v.corrected.size() != v.sequence.size()) {
       CHECK_LT(v.corrected.size(), v.sequence.size());
       v.sequence = v.sequence.subseq(0, v.corrected.size());
     }
     m_stats.corrected_read_count++;
     m_stats.corrected_read_bases += v.sequence.size();
     return true;
   }

   if (m_params.exact) {
     unsigned good_bases = verify_kmers(v.sequence, *m_kdb);
     unsigned needed_good_bases = m_params.trim_after_portion * v.sequence.size();
     if (good_bases >= needed_good_bases) {
       if (good_bases < v.sequence.size()) {
         v.sequence = v.sequence.subseq(0, good_bases);
       }
       m_stats.corrected_read_count++;
       m_stats.corrected_read_bases += v.sequence.size();
       // Double check to make sure no one uses length of
       // v.sequence; should use length of v.corrected instead.
       // Rename v.sequence to v.uncorrected maybe?
       v.corrected = v.sequence;
       return true;
     } else {
       m_stats.failed_correction_count++;
       return false;
     }
	}

	if (v.sequence.size() > m_stats.corrected_base_dist.size())
	{
		//LT_INFO
		// SPLOG("correct_reads_mapper::map_one_read> Resizing m_stats.corrected_base_dist to %lu", v.sequence.size());
		m_stats.corrected_base_dist.resize(v.sequence.size(), 0);
	}

	std::vector<kmer_t> path;
	double c = align_kmer(path, v.sequence, v.quality, *m_kdb, m_params.min_base_quality, m_params.max_quality_cost);
	if (c >= m_params.max_quality_cost)
	{
		m_stats.failed_correction_count++;
		if (m_params.trace)
			SPLOG("correct_reads_mapper::map_one_read> Unable to map: %s", id.pair_name.c_str());
		return false;
	}

	v.corrected = get_corrected(path, m_kdb->kmer_size());
	m_stats.corrected_read_count++;
	m_stats.corrected_read_bases += v.corrected.size();

	uint16_t base_diff_count = 0;
	auto corrected_iter = v.corrected.begin();
	auto original_read_iter = v.sequence.begin();
	while (corrected_iter != v.corrected.end())
	{
		if (*corrected_iter++ != *original_read_iter++) { base_diff_count++; }
	}

	m_stats.corrected_base_dist[base_diff_count]++;

	if (m_params.trace)
	{
		SPLOG("correct_reads_mapper::map_one_read> Read %s original sequence  = %s", id.pair_name.c_str(), v.sequence.as_string().c_str());
		SPLOG("correct_reads_mapper::map_one_read> Read %s corrected sequence = %s", id.pair_name.c_str(), v.corrected.as_string().c_str());
		v.trace_me = true;
	}
	return true;
}

void correct_reads_mapper::install_metadata1(meta::data& metadata)
{
	// SPLOG_P(LOG_DEBUG, "correct_reads_mapper::install_metadata>");
	truncate_corrected_base_dist(m_stats.corrected_base_dist);
	metadata.set(meta::ns::readonly, "corrected_read_count", m_stats.corrected_read_count);
	metadata.set(meta::ns::readonly, "corrected_read_bases", m_stats.corrected_read_bases);
	metadata.set(meta::ns::readonly, "corrected_base_dist", m_stats.corrected_base_dist);
	metadata.set(meta::ns::readonly, "failed_correction_count", m_stats.failed_correction_count);
}

void correct_reads_mapper::truncate_corrected_base_dist(std::vector<uint64_t>& dist_vector)
{
	// SPLOG_P(LOG_DEBUG, "correct_reads_mapper::truncate_corrected_base_dist>");
	auto last_non_zero_iter = std::find_if(dist_vector.rbegin(), dist_vector.rend(),
		[](uint64_t element) { return element != 0; }
	).base();
	if (last_non_zero_iter != dist_vector.end() && last_non_zero_iter != dist_vector.begin()) {
		dist_vector.erase(++last_non_zero_iter, dist_vector.end());
	}
}

#pragma once

#include "modules/bio_mapred/align_kmer.h"
#include "modules/mapred/dual_mapper.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_mapred/correct_reads.h"
#include "modules/bio_base/coverage_record.h"

class correct_reads_mapper: public typed_dual_mapper<correct_reads_mapper,
	read_id, unaligned_reads,
	std::string, corrected_reads,
	seq_position, coverage_record>
{
public:
	correct_reads_mapper(const std::string& params);

	void set_watchdog(const std::function<void()>& watchdog) override { m_watchdog = watchdog; }
	void setup() override;
	void typed_map(const read_id& key, const unaligned_reads& value);
	void install_metadata1(meta::data& metadata) override;
	task_requirements get_requirements() override;

private:
	bool map_one_read(corrected_read& v, const read_id& id, const unaligned_read& r);
	// Truncate trailing zeroes from distribution vector.
	void truncate_corrected_base_dist(std::vector<uint64_t>& dist_vector);

private:
	correct_reads_params m_params;
	read_correction_stats m_stats;
	std::function<void()> m_watchdog;
	std::unique_ptr<const kmer_set> m_kdb;
};

#pragma once

#include "modules/io/log.h"
#include "modules/io/transfer_object.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/task.h"

// NOTE: This is for the correct_reads step. For correct_reads_only, see correct_reads.h
struct read_correction_params
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(min_kmer_score, TF_STRICT);
		FIELD(min_base_quality, TF_STRICT);
		FIELD(max_quality_cost, TF_STRICT);
		FIELD(trim, TF_STRICT);
		FIELD(trace, false);
		FIELD(with_coverage, false);
		FIELD(skip_snps, false);
		FIELD(exact, false);
		FIELD(sys_err_thresh, float(0.0));
		FIELD(rnd_err_thresh, float(0.0));
		FIELD(overrep);
		FIELD(overrep);
		FIELD(trim_after_portion);
		FIELD(frc_max_corrections);
		FIELD(frc_min_good_run);
	}

	size_t min_kmer_score;
	double min_base_quality;
	double max_quality_cost;
	size_t trim = 0;
	bool trace = false;
	bool with_coverage = false;
	bool skip_snps = false;
	bool exact = false;
	float sys_err_thresh = 0.0;
	float rnd_err_thresh = 0.0;
	manifest overrep;
	double trim_after_portion = 1.0;
	unsigned frc_max_corrections = 2;
	unsigned frc_min_good_run = 2;

	void validate()
	{
		SPLOG_P(LOG_DEBUG,
			"read_correction_params::validate> min_kmer_score: %lu, min_base_quality: %0.2f, max_quality_cost: %0.2f, "
			"trim: %lu, trace: %s, skip_snps: %s, exact: %s",
			min_kmer_score,
			min_base_quality,
			max_quality_cost,
			trim,
			trace ? "true" : "false",
			skip_snps ? "true" : "false",
			exact ? "true" : "false"
		);
	}
};

class read_correction_task : public task_impl<read_correction_task>
{
public:
	static std::string s_type() { return "read_correction_task"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(reads, TF_STRICT);
		FIELD(kmers, TF_STRICT);
		FIELD(params, TF_STRICT);
		FIELD(m_state, TF_STRICT);
		FIELD(m_subtask, TF_STRICT);
		FIELD(kmers_filt, TF_STRICT);
		FIELD(kdb_man, TF_STRICT);
	}

	void run();

	manifest reads;
	manifest kmers;
	manifest kmers_filt;
	manifest kdb_man;
	read_correction_params params;

private:
	int m_state = 0; // 0 = nothing done, 1 = kmers filtered, 2 = kdb made, 3 = reads corrected
	subtask_id m_subtask;
};

#pragma once

#include "modules/mapred/manifest.h"
#include "modules/mapred/task.h"
#include "modules/io/transfer_object.h"
#include "modules/build_seqset/kmer_counter.h"

struct kmerize_bf_params
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(kmer_size, TF_STRICT);
		FIELD(partitions, TF_STRICT);
		FIELD(read_length, TF_STRICT);
		OBSOLETE_FIELD(trim, TF_STRICT);
		FIELD(read_parts, TF_STRICT);
		FIELD(overrep, size_t(1024*1024*1024));
		FIELD(prior_count, size_t(5));
		FIELD(skew_cutoff, float(0.0));
		FIELD(error_rate, TF_STRICT);
		FIELD(reference, TF_STRICT);
		FIELD(ref_size, TF_STRICT);
		FIELD(memory_bound);
		FIELD(num_threads);
		FIELD(min_count);
        FIELD(rnd_err_thresh);
        FIELD(sys_err_thresh);
        FIELD(dump_kmers_file);
	}

	size_t kmer_size = 0;
	size_t partitions = 0;
	size_t read_parts = 0;
	size_t read_length = 0;
	size_t overrep = 1024*1024*1024;
	size_t prior_count = 5;
	float skew_cutoff = 0.0;
	double error_rate;
	std::string reference;
	size_t ref_size = 0;
	size_t memory_bound = 0;
	size_t num_threads = 0;
	size_t min_count = 4;
  double rnd_err_thresh = 0;
  double sys_err_thresh = 0;
  std::string dump_kmers_file;

	void validate();
};

class kmerize_bf_subtask : public task_impl<kmerize_bf_subtask>
{
public:
	static std::string s_type() { return "kmerize_bf_subtask"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(params, TF_STRICT);
	}

	void run();

	task_requirements get_requirements() override
	{
		return task_requirements {
			.profile = "himem",
			.cpu_minutes = 60,
		};
	}

	manifest input;
	kmerize_bf_params params;
};

// Executes kmerize_bf_subtask, but uses an existing kmer counter.
// This kmer counter should already have executed the probabilistic
// pass.
std::pair<std::unique_ptr<kmer_set>, std::vector<manifest>> run_kmerize_subtask(
    const kmerize_bf_params& params, const manifest& input,
    build_seqset::kmer_counter* counter,
    progress_handler_t update_progress = null_progress_handler);
std::vector<std::string> get_kmer_filter_result_types();

class kmerize_bf_task : public task_impl<kmerize_bf_task>
{
public:
	static std::string s_type() { return "kmerize_bf_task"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(params, TF_STRICT);
		FIELD(m_state, TF_STRICT);
		FIELD(m_subtask, TF_STRICT);
		FIELD(m_histogram, TF_STRICT);
		FIELD(m_overrep, TF_STRICT);
		FIELD(m_kmer_counts);
	}

	void run();

	enum state
	{
		kmerize,
		sort_from_kmerize,
		final,

        // Start at the "sort" stage from externally provided kmer
        // counts instead of counting them internally.
        sort,

        // Does nothing, but supplies kmer_counts, histogram, and
        // overrep on output.
        do_nothing,
	};

	manifest input;
	kmerize_bf_params params;

	size_t m_state = state::kmerize;
	subtask_id m_subtask;
	manifest m_kmer_counts;
	manifest m_histogram;
	manifest m_overrep;
};

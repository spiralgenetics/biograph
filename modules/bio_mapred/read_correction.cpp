#include "modules/bio_mapred/read_correction.h"
#include "modules/bio_mapred/align_kmer.h"
#include "modules/bio_mapred/kmerize_reads_mapper.h"
#include "modules/bio_mapred/kmer_filter_mapper.h"
#include "modules/bio_mapred/kmerize_bf.h"
#include "modules/bio_mapred/correct_reads.h"
#include "modules/bio_mapred/kmers_to_db.h"
#include "modules/bio_mapred/compute_coverage.h"
#include "modules/bio_mapred/sv_call_reducer.h"
#include "modules/bio_mapred/read_pileup_reducer.h"
#include "modules/bio_base/reference.h"
#include "modules/mapred/map_reduce_task.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/map_task.h"
#include "modules/mapred/dual_map_task.h"
#include "modules/mapred/reduce_task.h"
#include "modules/mapred/sort_task.h"
#include "modules/io/make_unique.h"
#include "tools/version.h"

REGISTER_TASK(read_correction_task);

void read_correction_task::run() {
  SPLOG("read_correction_task::run> Start");
  if (m_state == 0 && kmers.metadata().has_key(meta::ns::internal, "kmer_db"))  {
    // We were passed in a kmer set; skip directly to read correction.
		split_progress(.01, .01);  // .01 for prep, .98 for children, .01 for future

		// Load KDB
		std::string kdb = kmers.metadata().get<std::string>(meta::ns::internal, "kmer_db");
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::2> Kmer DB record count = %lu",
            kdb_man.get_num_records());

		// Setup correct_reads_params
		correct_reads_params crp;
		crp.kmer_db = kdb;
		crp.min_base_quality = params.min_base_quality;
		crp.max_quality_cost = params.max_quality_cost;
		crp.skip_snps = params.skip_snps;
		crp.trim = params.trim;
		crp.exact = params.exact;
        crp.trim_after_portion = params.trim_after_portion;
        crp.frc_max_corrections = params.frc_max_corrections;
        crp.frc_min_good_run = params.frc_min_good_run;
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::2> trim %lu", params.trim);

		// Kick off correction
		auto task = make_unique<dual_map_task>();
		task->input = reads;
		task->map = "correct_reads";
		task->map_param = json_serialize(crp);
		m_subtask = add_subtask(std::move(task));
		m_state = 3;
  }
	else if (m_state == 0) {
		// Filter the kmers
		split_progress(.01, .95);  // .01 for prep, .04 for children, .95 for future
		kmer_filter_params kp;
		kp.min_count = params.min_kmer_score;
		kp.kmer_size = kmers.metadata().get<unsigned>(meta::ns::readonly, "kmer_size");
		kp.overrep = params.overrep;
		kp.sys_err_thresh = params.sys_err_thresh;
		kp.rnd_err_thresh = params.rnd_err_thresh;
		auto task = make_unique<map_task>();
		task->input = kmers;
		task->map = "kmer_filter";
		task->map_param = json_serialize(kp);
		task->stable_sort = true;
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::0> Prefiltered k-mer count: %lu", kmers.get_num_records());
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::0> k-mer filter params: \"%s\"", json_serialize(kp).c_str());
		m_subtask = add_subtask(std::move(task));
		m_state = 1;
	}
	else if (m_state == 1) {
		// Get results
		get_output(kmers_filt, m_subtask);
		kmers_filt.update_metadata(kmers);
		kmers_filt.metadata().set(meta::ns::readonly, "filtered_kmers", kmers.get_num_records() - kmers_filt.get_num_records());
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::1> Filtered k-mer record count = %zu", kmers_filt.get_num_records());
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::1> Filtered k-mer size = %zu", kmers_filt.get_num_records());
		split_progress(.01, .95);  // .01 for prep, .04 for children, .95 for future

		// Setup kmers_to_db
		auto task = make_unique<kmers_to_db_task>();
		task->input = kmers_filt;
		m_subtask = add_subtask(std::move(task));
		m_state = 2;
	}
	else if (m_state == 2) {
		split_progress(.01, .01);  // .01 for prep, .98 for children, .01 for future

		// Load KDB
		get_output(kdb_man, m_subtask);
		kdb_man.update_metadata(kmers_filt);
		std::string kdb = kdb_man.metadata().get<std::string>(meta::ns::internal, "kmer_db");
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::2> Kmer DB record count = %lu", kdb_man.get_num_records());

		// Setup correct_reads_params
		correct_reads_params crp;
		crp.kmer_db = kdb;
		crp.min_base_quality = params.min_base_quality;
		crp.max_quality_cost = params.max_quality_cost;
		crp.skip_snps = params.skip_snps;
		crp.trim = params.trim;
		crp.exact = params.exact;
        crp.trim_after_portion = params.trim_after_portion;
        crp.frc_max_corrections = params.frc_max_corrections;
        crp.frc_min_good_run = params.frc_min_good_run;
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::2> trim %lu", params.trim);

		// Kick off correction
		auto task = make_unique<dual_map_task>();
		task->input = reads;
		task->map = "correct_reads";
		task->map_param = json_serialize(crp);
		m_subtask = add_subtask(std::move(task));
		m_state = 3;
	}
	else if (m_state == 3) {
		// Load corrected reads
		std::vector<manifest> both;
		get_output(both, m_subtask);
		manifest& cr_man = both[0];
		cr_man.update_metadata(reads, kdb_man);
		SPLOG_P(LOG_DEBUG, "read_correction_task::run::3> Corrected reads count = %lu", cr_man.get_num_records());

		// Set as output
		if (params.with_coverage) {
			set_output(both);
		} else {
			set_output(cr_man);
		}
		SPLOG("read_correction_task::run> Done")
	}
}


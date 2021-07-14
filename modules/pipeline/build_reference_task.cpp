#include <boost/filesystem.hpp>

#include "modules/bio_base/bwt_file.h"
#include "modules/bio_base/flat_ref.h"
#include "modules/bio_format/fasta_ref_importer.h"
#include "modules/bio_mapred/make_bwt.h"
#include "modules/io/defaults.h"
#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/io/progress_tracker_types.h"
#include "modules/pipeline/build_reference_task.h"
#include "modules/pipeline/dataset_meta.h"
#include "modules/pipeline/dataset_path.h"
#include "tools/version.h"

REGISTER_TASK(build_reference_task);

namespace fs = boost::filesystem;

build_reference_task::build_reference_task(
	const std::string& out_dir,
	const std::string& reference_name)
	: m_out_dir(out_dir)
	, m_ref_name(reference_name)
{
	SPLOG("build_reference_task::build_reference_task> making reference for: %s", reference_name.c_str());
}

void build_reference_task::run()
{
	std::string prefix = "build_reference_task::run>";
	SPLOG("%s entering run, saving to %s", prefix.c_str(), m_out_dir.c_str());

	// TODO: implement real progress for flat and bwt generation.
	// These are estimates that do not smoothly increase!
	static const double import_fraction = 0.26;                 // 26% of the total time importing fasta
	static const double flat_fraction = import_fraction + 0.32; // 32% flat
	static const double bwt_fraction = flat_fraction + 0.39;    // 42% bwt, but underreport so we're not hung on 100%

	path out_dir_path(m_out_dir);
	path new_ref_dir_path = out_dir_path.append(m_ref_name);
	path fasta_tmp_path = new_ref_dir_path.append(defaults.original_fasta);

	path reference_fasta_path = new_ref_dir_path.append(defaults.reference_fasta);
	path alu_fasta_path = new_ref_dir_path.append(defaults.alu_fasta);
	path reference_ref_path = new_ref_dir_path.append(defaults.reference_ref);
	path reference_bwt_path = new_ref_dir_path.append(defaults.reference_bwt);

	size_t fasta_size = fasta_tmp_path.size();

	if (m_state == IMPORT_FASTA) {
		split_progress(import_fraction, 1 - import_fraction);
		// =========  IMPORT FASTA ===============================================
		SPLOG("%s importing new reference fasta file for %s", prefix.c_str(), m_ref_name.c_str());

		if (fasta_size == 0) {
			throw io_exception(prefix + " fasta file '" + m_ref_name + "' is empty!");
		}

		const size_t write_modulo = 1 + fasta_size / 20000;

		progress_t update = [&](size_t read, size_t written) {
			update_progress(written / (double)(2*fasta_size));
			return write_modulo;
		};

		std::vector<std::string> dummy_scaf_order;

		file_reader raw_in(fasta_tmp_path.bare_path());
		fasta_ref_importer fri(new_ref_dir_path.bare_path(), raw_in, dummy_scaf_order, m_min_n_run, update);
		fri.run();
		raw_in.close();

		update_progress(1.0);

		// =========  DONE IMPORTING FASTA ===============================================


		// ======== VERIFY FASTA IMPORT OUTPUT ==============================================
		std::string import_failed(prefix + " fasta import failed: ");
		path::exist_enum e = reference_fasta_path.exists();
		if (e != path::e_file) {
			throw io_exception(import_failed + reference_fasta_path.bare_path() + " does not exist.");
		}

		path ref_path = new_ref_dir_path.append(defaults.karyotype);
		if (ref_path.exists() == path::e_no_exist) {
			throw io_exception(import_failed + ref_path.bare_path() + " does not exist");
		}
		// ======== DONE VERIFYING FASTA IMPORT OUTPUT =======================================

		m_state = MAKE_FLAT;
	}
	else if (m_state == MAKE_FLAT) {
		SPLOG("import done");
		SPLOG("Making .ref");
		split_progress(flat_fraction, 1 - flat_fraction);

		// Make the flat reference from the pre-processed fasta.
		// This will change when the old reference_assembly stuff is removed.
		file_reader in(reference_fasta_path.bare_path());
		file_writer flat(reference_ref_path.bare_path());

		flat_ref_builder flb(in, flat);
		flb.run();

		try {
			flat_ref ref(reference_ref_path.bare_path());
		}
		catch(const io_exception&) {
			fs::remove(reference_ref_path.bare_path());
			throw;
		}
		m_state = MAKE_BWT;
		update_progress(1.0);
	}
	else if (m_state == MAKE_BWT) {
		SPLOG("Making .ref done");
		SPLOG("Making bwt");
		split_progress(bwt_fraction, 1 - bwt_fraction);
		auto bwt_task = make_unique<make_bwt_task>();
		bwt_task->input_ref = reference_ref_path.bare_path();
		bwt_task->output_bwt = reference_bwt_path.bare_path();

		add_subtask(std::move(bwt_task));

		update_progress(0.99);

		m_state = DONE;
	}
	else if (m_state == DONE) {
		SPLOG("Making bwt done");

		dataset_meta dm;
		dm.type = datatype_registry::find("reference");

		// successfully terminate the task by setting the output
		update_progress(1.0);
		set_output(dm.the_manifest);
	}
	else {
		throw io_exception(prefix + " Unknown state");
	}
}

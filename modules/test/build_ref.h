#pragma once

#include "modules/test/build_ref.h"
#include "modules/io/config.h"
#include "modules/io/defaults.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/pipeline/build_reference_task.h"
#include "modules/test/test_utils.h"
#include "tools/version.h"

#include <string>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;
struct build_ref_impl
{
	void setup()
	{
		std::string out_dir(CONF_S(reference_path));
		test_dir = fs::path(out_dir);
		task_path = test_dir / "task";

		fs::remove_all(test_dir);
	}

	void run_task(const std::string& ref_name, const std::string& fasta_path, const std::string& alu_fasta_path = "")
	{
		ASSERT_TRUE(fs::exists(fasta_path)) << fasta_path;
		fs::path ref_dir = test_dir / ref_name;
		fs::create_directories(ref_dir);
		fs::path fasta(fasta_path);
		// "fasta" may be a relative path into our golden directory,
		// so make it absolute before creating this symlink.
		fs::create_symlink(fs::canonical(fasta), ref_dir / defaults.original_fasta);

		if (! alu_fasta_path.empty()) {
			file_reader alu_fasta_reader{alu_fasta_path};
			file_writer alu_fasta_writer{(ref_dir / "alu.fasta").native()};
			io_copy(alu_fasta_reader, alu_fasta_writer);
			alu_fasta_writer.close();
		}

		tm.run_task(out, task_path.native(), make_unique<build_reference_task>(test_dir.native(), ref_name));
	}

	fs::path test_dir;
	fs::path task_path;
	manifest out;
	task_mgr_local tm;
};


void perform_build_ref(const std::string& ref_name, const std::string& fasta_path, const std::string& alu_fatsa_path = "");

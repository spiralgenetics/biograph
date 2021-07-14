#include "modules/test/build_ref.h"
#include "modules/io/config.h"
#include "modules/io/defaults.h"
#include "modules/io/file_io.h"
#include "modules/io/make_unique.h"
#include "modules/mapred/task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/pipeline/build_reference_task.h"
#include "modules/test/test_utils.h"
#include "tools/version.h"

#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

namespace fs = boost::filesystem;

void perform_build_ref(const std::string& ref_name, const std::string& fasta_path, const std::string& alu_fasta_path)
{
	build_ref_impl impl;
	impl.setup();
	impl.run_task(ref_name, fasta_path, alu_fasta_path);
}

#include "modules/test/test_utils.h"
#include "gtest/gtest.h"

#include "modules/io/file_io.h"
#include "modules/bio_format/fastq.h"
#include "modules/bio_mapred/phred64_to_33_mapper.h"
#include "modules/test/local_context.h"
#include "modules/test/fastq_test_utils.h" 

#include <boost/filesystem.hpp>

TEST(Phred, 64_to_33)
{
	boost::filesystem::create_directory(make_path("phred64_to_33"));
	
	local_context context(2, 1000000, path(make_path("phred64_to_33/local_context")));
	make_fastq_kv("golden/E_coli_phred64.fq", make_path("phred64_to_33/E_coli_phred64.kvp"));

	try {
		manifest phred64_manifest;
		phred64_manifest.add(file_info(path(make_path("phred64_to_33/E_coli_phred64.kvp")), 2290, 10), 0);

		manifest phred33_manifest = context.map_only("phred64_to_33", "", phred64_manifest);
		ASSERT_EQ(phred33_manifest.get_num_records(), static_cast<unsigned>(10));

		manifest_reader a_manifest_reader(phred33_manifest);
		kv_reader a_kv_reader(a_manifest_reader);
		auto out_path = make_path("phred64_to_33/E_coli_phred33.fq");
		path phred33_reads_path(out_path);
		std::unique_ptr<writable> phred33_writable_ptr(phred33_reads_path.write());

		fastq_exporter	the_fastq_exporter(*phred33_writable_ptr);
		the_fastq_exporter.export_from(a_kv_reader);

		ASSERT_TRUE(diff(out_path, "golden/E_coli_phred33.fq"));
	}
	catch (const io_exception& e)
	{
		printf("Phred, 64_to_33 Exception: %s\n", e.message().c_str());
		throw;
	}
}


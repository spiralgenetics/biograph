#include "modules/test/test_utils.h"
#include "modules/mapred/map_reduce_task.h"
#include "modules/bio_format/read_qual.h"
#include "modules/bio_format/fastq.h"
#include "modules/io/file_io.h"
#include "modules/io/zip.h"
#include "gtest/gtest.h"
#include "modules/mapred/task_mgr.h"
#include "modules/io/config.h"
#include "modules/mapred/sort_task.h"
#include "modules/bio_mapred/kmerize_reads_mapper.h"
#include "modules/bio_format/kmer_count.h"
#include "modules/mapred/ex_im_porter_data.h"
#include "modules/pipeline/paired_merger.h"
#include "modules/test/local_context.h"
#include "modules/test/fastq_test_utils.h"

TEST(pipeline, LittlePipeline)
{
	make_fastq_kv("golden/e_coli_10000snp.fq", make_path("e_coli_10000.kvp"));

	local_context context(1, 500000, make_path("little"));

	manifest e_coli_reads;
	e_coli_reads.add(file_info(path(make_path("e_coli_10000.kvp")), 1017780, 10000), 0);

	manifest qual = context.map_reduce("read_qual", "", "lexical", "sum", "", e_coli_reads, true);

	ASSERT_EQ(qual.get_num_records(), 141U);

	manifest_reader mr(qual);
	kv_reader kv(mr);
	path out_path(make_path("pipeline_test/simple.txt"));
	std::unique_ptr<writable> writer(out_path.write());

	read_qual_exporter exporter_(*writer);
	exporter_.export_from(kv);

	ASSERT_TRUE(diff(out_path, "golden/simple.txt"));
}

TEST(pipeline, PairedReadQual)
{
	ASSERT_NO_THROW(make_zipped_fastq_kv("golden/ftest/ERR_1.fastq.gz", make_path("ERR_1.kvp")));

	file_reader	ERR1_kvp_file_reader(make_path("ERR_1.kvp"));
	kv_reader	ERR1_kv_reader(ERR1_kvp_file_reader);

	file_reader ERR2_file_reader("golden/ftest/ERR_2.fastq.gz");
	zip_reader ERR2_unzipper(ERR2_file_reader);
	file_writer paired_file_writer(make_path("ERR_paired.kvp"));
	kv_writer paired_kv_writer(paired_file_writer);

	fastq_importer fastq_importer(ERR2_unzipper);
	paired_merger	ERR_merge_writer(paired_kv_writer, ERR1_kv_reader);
	fastq_importer.import(ERR_merge_writer, discard_simple_metadata());
	ERR2_file_reader.close();
	paired_file_writer.close();
	ERR1_kvp_file_reader.close();

	local_context context(1, 500000, make_path("read_qual"));

	manifest ERR1_paired;
	ERR1_paired.add(file_info(path(make_path("ERR_paired.kvp")), 17957421, 60000), 0);

	manifest qual = context.map_reduce("read_qual", "", "lexical", "sum", "", ERR1_paired, true);

	manifest_reader qual_manifest_reader(qual);
	kv_reader read_qual_kv_reader(qual_manifest_reader);
	path out_path(make_path("read_qual_test/ERR1_read_qual.json"));
	std::unique_ptr<writable> writer(out_path.write());

	read_qual_exporter exporter_(*writer);
	exporter_.export_from(read_qual_kv_reader);

	ASSERT_TRUE(diff(out_path, "golden/ERR_read_qual.json"));
}

TEST(pipeline, KmerPipeline)
{
	make_fastq_kv("golden/e_coli_10000snp.fq", make_path("e_coli_10000.kvp"));

	local_context context(1, 500000, make_path("kmer"));

	manifest e_coli_reads;
	e_coli_reads.add(file_info(path(make_path("e_coli_10000.kvp")), 1017780, 10000), 0);

	kmerize_reads_params kp;
	std::string params = R"|(
		{
			"kmer_size" : 23,
			"trim" : 0,
			"use_score" : false
		}
	)|";
	json_deserialize(kp, params);
	kp.validate();

	manifest kmer = context.map_reduce("kmerize_reads",
		json_serialize(kp), "lexical", "kcount", "", e_coli_reads, true);

	ASSERT_EQ(kmer.get_num_records(), 9978U);

	manifest_reader mr(kmer);
	kv_reader kv(mr);
	path out_path(make_path("kmers.txt"));
	std::unique_ptr<writable> writer(out_path.write());

	kmer_count_exporter exporter_(*writer, kp.kmer_size);
	exporter_.export_from(kv);

	ASSERT_TRUE(diff(out_path, "golden/kmers.txt"));
}

TEST(pipeline, kmer_count)
{
	const auto reads_kvp = make_path("quick_e_coli.kvp");
	make_fastq_kv("golden/quick_e_coli.fq", reads_kvp);

	local_context context(1, 500000, make_path("kmer_count"));

	manifest reads;
	reads.add(file_info(path(reads_kvp), 480171, 2223), 0);

	kmerize_reads_params kp;
	std::string params = R"|(
		{
			"kmer_size" : 23,
			"trim" : 0,
			"use_score" : false
		}
	)|";
	json_deserialize(kp, params);
	kp.validate();

	manifest kmers = context.map_reduce("kmerize_reads",
		json_serialize(kp), "lexical", "kcount", "", reads, true);

	ASSERT_EQ(kmers.get_num_records(), 29397U);

	manifest_reader mr(kmers);
	kv_reader kvr(mr);
	path out_path(make_path("kmer_count.txt"));
	std::unique_ptr<writable> writer(out_path.write());
	kmer_count_exporter exporter_(*writer, kp.kmer_size);
	exporter_.export_from(kvr);

	ASSERT_TRUE(diff(out_path, "golden/kmer_count.txt"));
}

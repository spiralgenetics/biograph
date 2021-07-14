#include "modules/test/test_utils.h"
#include "modules/mapred/task_mgr.h"
#include "modules/test/fastq_test_utils.h"
#include "modules/bio_mapred/kmerize_reads_mapper.h"
#include "modules/mapred/map_reduce_task.h"
#include <gtest/gtest.h>


TEST(read_trim, basic)
{
	task_mgr_local tm;
	path out_path(make_path("read_trim_test"));

	// Get ecoli data
	manifest e_coli_reads;
	make_fastq_kv("golden/e_coli_10000snp.fq", make_path("e_coli_10000.kvp"));
	e_coli_reads.add(file_info(path(make_path("e_coli_10000.kvp")), 1017780, 10000), 0);

	// Kmerize reads
	// No-trim run
	kmerize_reads_params kp_no_trim;
	json_deserialize(kp_no_trim, R"|(
			{
				"kmer_size" : 30,
				"trim" : 0,
				"use_score" : false
			}
		)|"
	);

	std::unique_ptr<map_reduce_task> t = make_unique<map_reduce_task>();
	t->input = e_coli_reads;
	t->map = "kmerize_reads";
	t->map_param = json_serialize(kp_no_trim);
	t->sort = "lexical";
	t->reduce = "kcount";
	t->is_summary = true;
	t->use_sort = true;
	manifest kmers;
	tm.run_task(kmers, out_path.append("kmers"), std::move(t));
	auto no_trim_records = kmers.get_num_records();

	// Trim run
	kmerize_reads_params kp_trim;
	json_deserialize(kp_trim, R"|(
			{
				"kmer_size" : 30,
				"trim" : 32,
				"use_score" : false
			}
		)|"
	);

	std::unique_ptr<map_reduce_task> s = make_unique<map_reduce_task>();
	s->input = e_coli_reads;
	s->map = "kmerize_reads";
	s->map_param = json_serialize(kp_no_trim);
	s->sort = "lexical";
	s->reduce = "kcount";
	s->is_summary = true;
	s->use_sort = true;
	s->map_param = json_serialize(kp_trim);
	tm.run_task(kmers, out_path.append("kmers"), std::move(s));
	auto trim_records = kmers.get_num_records();

	ASSERT_GT(no_trim_records, trim_records);
}


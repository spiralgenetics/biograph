#include "modules/test/test_utils.h"
#include "modules/bio_mapred/align_kmer.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/mapred/task_mgr.h"
#include "modules/test/fastq_test_utils.h"
#include "modules/bio_mapred/kmerize_reads_mapper.h"
#include "modules/mapred/map_reduce_task.h"
#include "modules/bio_mapred/kmer_filter_mapper.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"

#include <gtest/gtest.h>
#include <chrono>

TEST(align_kmer, test)
{
	task_mgr_local tm;
	path out_path(make_path("align_kmer_test"));

	// Get ecoli data
	manifest e_coli_reads;
	make_fastq_kv("golden/e_coli_10000snp.fq", make_path("e_coli_10000.kvp"));
	e_coli_reads.add(file_info(path(make_path("e_coli_10000.kvp")), 1017780, 10000), 0);

	// Kmerize reads
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

	std::unique_ptr<map_reduce_task> t = make_unique<map_reduce_task>();
	t->input = e_coli_reads;
	t->map = "kmerize_reads";
	t->map_param = json_serialize(kp);
	t->sort = "lexical";
	t->reduce = "kcount";
	t->is_summary = true;
	t->use_sort = true;
	manifest kmers;
	tm.run_task(kmers, out_path.append("kmers"), std::move(t));

	// Filter kmer reads
	kmer_filter_params kfp;
	params = R"|(
		{
			"min_count" : 3,
			"kmer_size" : 23
		}
	)|";
	json_deserialize(kfp, params);
	kfp.validate();

	std::unique_ptr<map_task> t2 = make_unique<map_task>();
	t2->input = kmers;
	t2->map = "kmer_filter";
	t2->map_param = json_serialize(kfp);
	t2->stable_sort = true;
	manifest filtered_kmers;
	tm.run_task(filtered_kmers, out_path.append("filtered_kmers"), std::move(t2));
	SPLOG("Filtered kmers = %d", int(filtered_kmers.get_num_records()));

	// TODO: Clean up location of test data
	manifest_reader mr(filtered_kmers);
	kmer_set ks(mr, filtered_kmers.get_num_records(), kp.kmer_size);
	dna_sequence read("CCGGCGGTGACACCTGTTGATGGTGCATAGCTCGG");
	std::string  qual("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
	std::vector<kmer_t> out;
	double r = align_kmer(out, read, qual, ks, 20.0, 100);
	printf("r = %f\n", r);
	ASSERT_EQ(72.0, r);
	ASSERT_GT(out.size(), 1UL);
	for (size_t i = 0; i < out.size(); i++) {
		printf("%s\n", dna_sequence(out[i], ks.kmer_size()).as_string().c_str());
	}
	dna_sequence corrected = get_corrected(out, ks.kmer_size());
	printf("Orig = %s\n", read.as_string().c_str());
	printf("Corr = %s\n", corrected.as_string().c_str());
	dna_sequence read_fixed("GCGGCGGTGACACCTGTTGATGGTGCATTGCTCGG");
	ASSERT_EQ(corrected.as_string(), read_fixed.as_string());
}

// This test is a copy of the basic align_kmer test with a larger dataset for use in
// benchmarking the kmerization.  Normally the test is disabled.
TEST(align_kmer, DISABLED_benchmark)
{
	task_mgr_local tm;
	path out_path(make_path("align_kmer_benchmark"));

	manifest reads;
	make_zipped_fastq_kv("golden/ftest/ERR_1.fastq.gz", make_path("ERR_1.kvp"));
	reads.add(file_info(path(make_path("ERR_1.kvp")), 9191614, 30000), 0);

	// Kmerize reads
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

	std::unique_ptr<map_reduce_task> t = make_unique<map_reduce_task>();
	t->input = reads;
	t->map = "kmerize_reads";
	t->map_param = json_serialize(kp);
	t->sort = "lexical";
	t->reduce = "kcount";
	t->is_summary = true;
	t->use_sort = true;
	manifest kmers;

	auto start = std::chrono::system_clock::now();
	tm.run_task(kmers, out_path.append("kmers"), std::move(t));
	std::chrono::duration<double> kmerization_time = std::chrono::system_clock::now() - start;
	std::cout << "Kmerization ran in " << kmerization_time.count() << " seconds.\n";
}


TEST(align_kmer, table)
{
	for (unsigned i = 0; i < 127; i++)
	{
		std::cout << i << ' ' << kmerize_reads_mapper::mg_log_lookup_table[i] << '\n';
	}
}


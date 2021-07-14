#include "modules/bio_mapred/kmerize_bf.h"
#include "modules/bio_base/kmer.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_format/fastq.h"
#include "modules/bio_format/kmer_count.h"
#include "modules/mapred/histogram_export.h"
#include "modules/mapred/task_mgr.h"
#include "modules/io/hash.h"
#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/io/bloom_filter.h"
#include "modules/io/zip.h"
#include "modules/io/stopwatch.h"

#include "modules/test/test_utils.h"
#include "modules/test/fastq_test_utils.h"
#include "gtest/gtest.h"

#include <fstream>

struct bloom_stats
{
	bloom_stats()
		: true_positives(0)
		, true_negatives(0)
		, false_positives(0)
		, false_negatives(0)
	{}

	void score(bool actual, bool expected)
	{
		if (expected) {
			if (actual) {
				true_positives++;
			} else {
				false_negatives++;
			}
		} else {
			if (actual) {
				false_positives++;
			} else {
				true_negatives++;
			}
		}
	}

	void analyze(double error_rate)
	{
		double false_positive_rate = (double)false_positives /
			(false_positives + true_negatives);

		SPLOG("True positives:      %7zu", true_positives);
		SPLOG("True negatives:      %7zu", true_negatives);
		SPLOG("False positives:     %7zu", false_positives);
		SPLOG("False negatives:     %7zu", false_negatives);
		SPLOG("False positive rate: %.4f", false_positive_rate);

		ASSERT_EQ(0, false_negatives);
		ASSERT_LT(false_positive_rate, error_rate);
	}

	size_t true_positives;
	size_t true_negatives;
	size_t false_positives;
	size_t false_negatives;
};

TEST(bloom, basic)
{
	const size_t CAPACITY = 100;
	const double ERROR_RATE = 0.01;

	bloom_filter<2> bf(CAPACITY, ERROR_RATE);
	for (size_t i = 0; i < (CAPACITY - CAPACITY * 0.10); i++) {
		bf.add(prime_hasher(i));
	}

	bloom_stats stats;
	for (size_t i = 0; i < CAPACITY * 2; i++) {
		auto count = bf.lookup(prime_hasher(i));
		// SPLOG("[%zu]: %zu", i, count);
		stats.score(count, i < (CAPACITY - CAPACITY * 0.10));
	}

	stats.analyze(ERROR_RATE);
}

template <typename Hasher, size_t KmerSize>
void test_kmer(const std::string& input_path)
{
	const size_t KMER_SIZE = KmerSize;
	const size_t CAPACITY = 2423184;
	const double ERROR_RATE = 0.05;

	bloom_filter<2> bf(CAPACITY, ERROR_RATE);
	SPLOG("cells  (m): %zu", bf.cells());
	SPLOG("hashes (k): %zu", bf.hashes());
	SPLOG("bitmap memory consumed: %zu bytes", bf.bitmap().memory_used());

	file_reader fin(input_path);
	zip_reader unzip(fin);
	fastq_reader reader(unzip);

	std::map<kmer_t, size_t> kmer_count;
	size_t counter = 0;

	auto delta = stopwatch([&] {
		read_id key;
		unaligned_reads value;
		while (reader.read_msgpack(key, value)) {
			for (const auto& read : value) {
				for (auto kmer : kmer_str_view(read.sequence, KMER_SIZE)) {
					bf.add(Hasher(kmer));
					counter++;
					kmer_count[kmer]++;
				}
			}
		}
	});

	SPLOG("Processed %zu kmers in %ld ms", counter, delta.count());
	SPLOG("Unique kmers: %zu", kmer_count.size());
	SPLOG("kmer_count memory usage: ~%zu bytes", kmer_count.size() * (sizeof(kmer_t) + sizeof(size_t)));

	size_t expected[4] = {0, 0, 0, 0};
	size_t actual[4] = {0, 0, 0, 0};

	std::map<kmer_t, size_t> bloom_count;
	for (const auto& item : kmer_count) {
		bloom_count[item.first] = bf.lookup(Hasher(item.first));
	}

	size_t expected_total = 0;
	for (const auto& item : kmer_count) {
		expected_total += item.second;
		if (item.second < 3) {
			expected[item.second]++;
		}
		else {
			expected[3]++;
		}
	}

	size_t actual_total = 0;
	for (const auto& item : bloom_count) {
		actual_total += item.second;
		actual[item.second]++;
	}

	for (size_t i = 0; i < 4; i++) {
		double error = 0.0;
		if (expected[i]) {
			double diff = (double)actual[i] - (double)expected[i];
			error = diff / (double)expected[i];
		}
		SPLOG("kmers[%zu] expected: %10zu actual: %10zu error: %.4f%%",
			i, expected[i], actual[i], error * 100.0);
	}
	double diff = (double)actual_total - (double)expected_total;
	double error_total = diff / (double)expected_total;
	SPLOG("total:   expected: %10zu actual: %10zu error: %.4f%%",
		expected_total, actual_total, error_total * 100.0);

	ASSERT_LT(error_total, ERROR_RATE);
}

TEST(bloom, kmer_prime)
{
	test_kmer<prime_hasher, 30>("golden/ftest/human_reads.fastq.gz");
}

TEST(bloom, task)
{
	const auto FILE_SIZE = 480171;
	const auto NUM_RECORDS = 2223;
	const auto reads_kvp = make_path("quick_e_coli.kvp");
	make_fastq_kv("golden/quick_e_coli.fq", reads_kvp);

	manifest reads;
	reads.add(file_info(path(reads_kvp), FILE_SIZE, NUM_RECORDS), 0);

	auto task = make_unique<kmerize_bf_task>();
	task->input = reads;

	json_deserialize(task->params, R"|(
			{
				"kmer_size" : 23,
				"ref_size" : 10240,
				"partitions" : 0,
				"read_length" : 0,
				"trim" : 0,
				"read_parts" : 0,
				"error_rate" : 0.10,
				"reference" : ""
			}
		)|"
	);
	task->params.validate();

	path task_path(make_path("bloom_task"));
	task_mgr_local tmgr;
	std::vector<manifest> out;
	tmgr.run_task(out, task_path, std::move(task));

	manifest_reader mr1(out[0]);
	kv_reader kvr1(mr1);
	path out_path1(task_path.append("kmer_count.txt"));
	std::unique_ptr<writable> writer1(out_path1.write());
	kmer_count_exporter count_exporter(*writer1, 23);
	count_exporter.export_from(kvr1);

	manifest_reader mr2(out[1]);
	kv_reader kvr2(mr2);
	path out_path2(task_path.append("histogram.txt"));
	std::unique_ptr<writable> writer2(out_path2.write());
	histogram_exporter histogram_exporter(*writer2);
	histogram_exporter.export_from(kvr2);

	// TODO(nils): Do we still need to produce a kmer count output somewhere?
	// ASSERT_TRUE(diff(out_path1, "golden/kmer_count_bf.txt"));
}

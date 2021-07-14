#include "gtest/gtest.h"

#include "modules/io/log.h"
#include "modules/io/file_io.h"
#include "modules/mapred/temp_file.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_format/dna_io.h"
#include "modules/bio_mapred/flatten_seqset.h"
#include "modules/bio_mapred/merge_flat_seqset.h"
#include "modules/test/test_utils.h"

TEST(seqset_merge, partition)
{
	std::vector<std::string> seqsets {"datasets/hiv/biograph/ERR732132.bg/seqset"};

	int num_threads = 4;
	flatten_seqset flattener(seqsets, num_threads);
	for (int i = 0; i < num_threads; i++) {
		ASSERT_EQ(flattener.find_partition_sequence(i, 4).as_string()[0], char(dna_base(i)));
		ASSERT_EQ(flattener.find_partition_sequence(i, 4).size(), 1);
	}

	flatten_seqset flattener8(seqsets, 8);
	ASSERT_EQ(flattener8.find_partition_sequence(0, 8).as_string(), "AA");
	ASSERT_EQ(flattener8.find_partition_sequence(0, 8).size(), 2);

	ASSERT_EQ(flattener8.find_partition_sequence(3, 8).as_string(), "CG");
	ASSERT_EQ(flattener8.find_partition_sequence(3, 8).size(), 2);

	ASSERT_EQ(flattener8.find_partition_sequence(7, 8).as_string(), "TG");
	ASSERT_EQ(flattener8.find_partition_sequence(7, 8).size(), 2);

	flatten_seqset flattener16(seqsets, 16);
	ASSERT_EQ(flattener16.find_partition_sequence(0, 16).as_string(), "AA");
	ASSERT_EQ(flattener16.find_partition_sequence(0, 16).size(), 2);

	ASSERT_EQ(flattener16.find_partition_sequence(6, 16).as_string(), "CG");
	ASSERT_EQ(flattener16.find_partition_sequence(6, 16).size(), 2);

	ASSERT_EQ(flattener16.find_partition_sequence(15, 16).as_string(), "TT");
	ASSERT_EQ(flattener16.find_partition_sequence(15, 16).size(), 2);

	flatten_seqset flattener32(seqsets, 32);
	ASSERT_EQ(flattener32.find_partition_sequence(0, 32).as_string(), "AAA");
	ASSERT_EQ(flattener32.find_partition_sequence(0, 32).size(), 3);

	ASSERT_EQ(flattener32.find_partition_sequence(12, 32).as_string(), "CGA");
	ASSERT_EQ(flattener32.find_partition_sequence(12, 32).size(), 3);

	ASSERT_EQ(flattener32.find_partition_sequence(31, 32).as_string(), "TTG");
	ASSERT_EQ(flattener32.find_partition_sequence(31, 32).size(), 3);

	flatten_seqset flattener64(seqsets, 64);
	ASSERT_EQ(flattener64.find_partition_sequence(0, 64).as_string(), "AAA");
	ASSERT_EQ(flattener64.find_partition_sequence(0, 64).size(), 3);

	ASSERT_EQ(flattener64.find_partition_sequence(24, 64).as_string(), "CGA");
	ASSERT_EQ(flattener64.find_partition_sequence(24, 64).size(), 3);

	ASSERT_EQ(flattener64.find_partition_sequence(63, 64).as_string(), "TTT");
	ASSERT_EQ(flattener64.find_partition_sequence(63, 64).size(), 3);
}

TEST(seqset_merge, power2)
{
	ASSERT_TRUE(flatten_seqset::is_power_of_2(1));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(2));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(4));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(8));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(16));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(32));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(64));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(128));
	ASSERT_TRUE(flatten_seqset::is_power_of_2(256));

	ASSERT_FALSE(flatten_seqset::is_power_of_2(5));
	ASSERT_FALSE(flatten_seqset::is_power_of_2(344));
	ASSERT_FALSE(flatten_seqset::is_power_of_2(44));

	ASSERT_FALSE(flatten_seqset::is_power_of_2(0));

	ASSERT_FALSE(flatten_seqset::is_power_of_2(-1));
	ASSERT_FALSE(flatten_seqset::is_power_of_2(-2));
	ASSERT_FALSE(flatten_seqset::is_power_of_2(-3));
	ASSERT_FALSE(flatten_seqset::is_power_of_2(-4));
}

TEST(seqset_merge, power4)
{
	ASSERT_TRUE(flatten_seqset::is_power_of_4(1));
	ASSERT_TRUE(flatten_seqset::is_power_of_4(4));
	ASSERT_TRUE(flatten_seqset::is_power_of_4(16));
	ASSERT_TRUE(flatten_seqset::is_power_of_4(64));
	ASSERT_TRUE(flatten_seqset::is_power_of_4(256));

	ASSERT_FALSE(flatten_seqset::is_power_of_4(5));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(344));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(44));

	ASSERT_FALSE(flatten_seqset::is_power_of_4(0));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(2));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(8));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(32));

	ASSERT_FALSE(flatten_seqset::is_power_of_4(-1));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(-2));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(-3));
	ASSERT_FALSE(flatten_seqset::is_power_of_4(-4));
}

// Merge the SEQSET files in the second argument and write the sequences as text to a single
// file named by the third argument.
void merge_seqsets_test(int num_threads, const std::vector<std::string>& seqset_paths, const std::string& merged_seqs_output_path)
{
	flatten_seqset flattener(seqset_paths, num_threads);
	std::multimap<int, std::shared_ptr<scoped_temp_file>> temp_file_map = flattener();
	SPLOG("Multimap size %lu", temp_file_map.size());
	ASSERT_EQ(temp_file_map.size(), seqset_paths.size() * num_threads);

	for (auto map_value : temp_file_map) {
		dna_reader a_dna_reader(make_unique<file_reader>(map_value.second->path()));
		dna_sequence partition_seq(flattener.find_partition_sequence(map_value.first, num_threads));
		dna_sequence the_sequence;
		for (int i = 0; i < 1000; i++) {
			the_sequence = a_dna_reader.read();
			if (the_sequence.size() == 0) break;
			if (flatten_seqset::is_power_of_4(num_threads)) {
				ASSERT_EQ(partition_seq.as_string(), the_sequence.subseq(0, partition_seq.size()).as_string());
			} else {
				ASSERT_EQ(partition_seq.subseq(0, partition_seq.size() - 1).as_string()
					, the_sequence.subseq(0, partition_seq.size() - 1).as_string());
				ASSERT_TRUE(partition_seq[partition_seq.size() - 1] == the_sequence[partition_seq.size() - 1]
					|| int(partition_seq[partition_seq.size() - 1]) == int(the_sequence[partition_seq.size() - 1]) - 1
				);
			}
		}
	}

	SPLOG("About to merge...");
	merge_flat_seqsets merger;
	std::vector<std::shared_ptr<scoped_temp_file>> merged_files = merger.merge_seqs(temp_file_map);
	SPLOG("Done merging, %lu files returned.", merged_files.size());

	file_writer merged_strings(merged_seqs_output_path);
	for (const auto temp_file : merged_files) {
		dna_reader temp_dna_reader(make_unique<file_reader>(temp_file->path()));
		dna_sequence the_seq = temp_dna_reader.read();
		 while (the_seq.size()) {
			merged_strings.write(the_seq.as_string().data(), the_seq.size());
			merged_strings.write("\n", 1);
			the_seq = temp_dna_reader.read();
		}
	}
	merged_strings.close();
}

TEST(seqset_merge, merge_one_text)
{
	std::vector<std::string> seqsets {"datasets/hiv/biograph/ERR381524.bg/seqset"};
	std::string merged_file_path {make_path("ERR381524.merge.seq")};
	merge_seqsets_test(32, seqsets, merged_file_path);

	std::string seqset_seq_out_path { make_path("ERR381524.seqset.seq") };
	int sys_return = ::system((boost::format("sh -c 'modules/biograph/bgbinary query --in %1% --query \"\" --verbose > %2%'")
		% seqsets[0] % seqset_seq_out_path).str().c_str());
	int wait_return = WEXITSTATUS(sys_return);
	ASSERT_EQ(wait_return, 0);
	ASSERT_TRUE(diff(seqset_seq_out_path, merged_file_path));
}

TEST(seqset_merge, threads)
{
	std::vector<std::string> seqsets {"datasets/hiv/biograph/ERR381524.bg/seqset"
		, "datasets/hiv/biograph/ERR732129.bg/seqset"
		, "datasets/hiv/biograph/ERR732130.bg/seqset"
		, "datasets/hiv/biograph/ERR732131.bg/seqset"
		, "datasets/hiv/biograph/ERR732132.bg/seqset"
	};

	std::vector<int> thread_counts = { 4, 8, 16, 32 };
	auto first_thread_count = thread_counts.front();

	std::string merged_file_path { make_path("merged_seqs") };
	for (const auto thread_num : thread_counts) {
		SPLOG("Merging with %d threads...", thread_num);
		merge_seqsets_test(thread_num, seqsets, merged_file_path + std::to_string(thread_num));

		if (thread_num != first_thread_count) {
			ASSERT_TRUE(diff(merged_file_path + std::to_string(first_thread_count)
				, merged_file_path + std::to_string(thread_num)));
		}
	}
}

std::vector<std::shared_ptr<scoped_temp_file>> make_temp_files(
	const std::vector<std::string>& seqset_paths
	, int num_threads
)
{
	SPLOG("Flattening SEQSETs...");
	flatten_seqset flattener(seqset_paths, 32);
	std::multimap<int, std::shared_ptr<scoped_temp_file>> temp_file_map = flattener();

	SPLOG("Merging flat files...");
	merge_flat_seqsets merger;
	return merger.merge_seqs(temp_file_map);
}

TEST(seqset_merge, multi_file_walk)
{
	std::vector<std::string> seqset_paths {"datasets/hiv/biograph/ERR381524.bg/seqset"
		, "datasets/hiv/biograph/ERR732129.bg/seqset"
		, "datasets/hiv/biograph/ERR732130.bg/seqset"
		, "datasets/hiv/biograph/ERR732131.bg/seqset"
		, "datasets/hiv/biograph/ERR732132.bg/seqset"
	};

	std::string merged_seqs_path = {make_path("merged_seqs")};
	SPLOG("Merging SEQSETs to a single file.");
	merge_seqsets_test(32, seqset_paths, merged_seqs_path);

	std::vector<std::shared_ptr<scoped_temp_file>> merged_files{make_temp_files(seqset_paths, 32)};

	SPLOG("Starting walk...");
	std::string walk_results_path = {make_path("walk_results")};
	file_writer walk_results(walk_results_path);
	multi_file_dna_buffer walker(std::move(merged_files));
	while(! walker.at_eof()) {
		walk_results.write(walker.get_sequence().as_string().data(), walker.get_sequence().size());
		walk_results.write("\n", 1);
		walker.advance();
	}
	walk_results.close();

	SPLOG("Diff of %s and %s", merged_seqs_path.c_str(), walk_results_path.c_str());
	ASSERT_TRUE(diff(merged_seqs_path, walk_results_path));
}

TEST(seqset_merge, merge_one)
{
	std::vector<std::string> seqsets {"datasets/hiv/biograph/ERR381524.bg/seqset"};
	std::string merged_file_path {make_path("ERR381524.gbwt")};

	SPLOG("Flattening SEQSET...");
	flatten_seqset flattener(seqsets, 32);
	std::multimap<int, std::shared_ptr<scoped_temp_file>> temp_file_map = flattener();

	SPLOG("Merging SEQSET...");
	merge_flat_seqsets()(merged_file_path, temp_file_map, true, 255 /* max read len */);

	SPLOG("Dumping sequences...");
	std::string seqset_seq_out_path { make_path("ERR381524.merged.seqset.seq") };
	std::string original_seq_out_path { make_path("ERR381524.original.seqset.seq") };
	int sys_return = ::system((boost::format("sh -c 'modules/biograph/bgbinary query --in %1% --query \"\" --verbose > %2%'")
		% seqsets[0] % original_seq_out_path).str().c_str());
	int wait_return = WEXITSTATUS(sys_return);
	ASSERT_EQ(wait_return, 0);

	sys_return = ::system((boost::format("sh -c 'modules/biograph/bgbinary query --in %1% --query \"\" --verbose > %2%'")
		% merged_file_path % seqset_seq_out_path).str().c_str());
	wait_return = WEXITSTATUS(sys_return);
	ASSERT_EQ(wait_return, 0);

	SPLOG("Comparing sequences...");
	ASSERT_TRUE(diff(original_seq_out_path, seqset_seq_out_path));
}

TEST(seqset_merge, merge_five)
{
	std::vector<std::string> seqsets {"datasets/hiv/biograph/ERR381524.bg/seqset"
		, "datasets/hiv/biograph/ERR732129.bg/seqset"
		, "datasets/hiv/biograph/ERR732130.bg/seqset"
		, "datasets/hiv/biograph/ERR732131.bg/seqset"
		, "datasets/hiv/biograph/ERR732132.bg/seqset"
	};
	std::string merged_file_path {make_path("MergedHIV.gbwt")};

	SPLOG("Flattening SEQSET...");
	flatten_seqset flattener(seqsets, 32);
	std::multimap<int, std::shared_ptr<scoped_temp_file>> temp_file_map = flattener();

	SPLOG("Merging SEQSET...");
	merge_flat_seqsets()(merged_file_path, temp_file_map, true, 255 /* max read len */);

	seqset_file merged_seqset_file(merged_file_path);
}

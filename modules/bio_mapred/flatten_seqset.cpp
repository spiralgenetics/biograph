#include <atomic>
#include <iostream>
#include <thread>

#include <boost/format.hpp>

#include "modules/bio_base/seqset.h"
#include "modules/bio_format/dna_io.h"
#include "modules/bio_mapred/flatten_seqset.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/utils.h"
#include "modules/mapred/temp_file.h"

flatten_seqset::flatten_seqset(std::vector<std::string> seqset_files, int num_threads, int max_read_size)
	: m_seqset_files(std::move(seqset_files))
	, m_num_threads(num_threads)
	, m_max_read_size(max_read_size)
{
	if (num_threads < 4) {
		throw io_exception(boost::format(
			"flatten_seqset requires a minimum of 4 threads, you"
			" requested %1%") % num_threads);
	}
	if (! is_power_of_2(num_threads)) {
		throw io_exception(boost::format(
			"flatten_seqset requires a number of threads that is a power of two, you"
			" requested %1%") % num_threads);
	}
}

std::multimap<int, std::shared_ptr<scoped_temp_file>> flatten_seqset::operator() ()
{
	SPLOG("flatten_seqset::operator()> Beginning flattening...");
	std::multimap<int, std::shared_ptr<scoped_temp_file>> temp_files_map;

	std::vector<std::future<void>> are_threads_done;

	for (const auto& seqset_path : m_seqset_files) {
		SPLOG("Opening seqset \"%s\".", seqset_path.c_str());
		std::cout << seqset_path.c_str() << std::endl;
		print_progress(0.0);

		seqset_file the_seqset_file(seqset_path);
		const seqset& the_seqset = the_seqset_file.get_seqset();
		SPLOG("Main: First shared = %u, context = %u", the_seqset.entry_shared(0), the_seqset.entry_size(0));
		the_seqset.populate_pop_front_cache();
		SPLOG("seqset \"%s\" opened and cache populated.", seqset_path.c_str());

		m_cur_progress = 0;
		size_t total_progress = the_seqset.size();

		for (int i = 0; i < m_num_threads; i++) {
			are_threads_done.emplace_back(
				std::async(std::launch::async,
					[this, i, &temp_files_map, &the_seqset]() { this->flatten_partition(i, temp_files_map, the_seqset); }
				)
			);
		}

		// Print progress in a separate thread.
		auto progress = std::async(std::launch::async,
			[&](){
				while(m_cur_progress < total_progress) {
					print_progress(float(m_cur_progress) / float(total_progress));
					usleep(400000);
				}
			}
		);

		for (int i = 0; i < m_num_threads; i++) {
			are_threads_done[i].get();
		}

		// Exit the progress thread
		m_cur_progress = total_progress;
		progress.get();

		are_threads_done.clear();
		print_progress(1.0);
		std::cout << std::endl;
	}

	SPLOG("flatten_seqset::operator()> Done flattening.");
	return temp_files_map;
}

void flatten_seqset::flatten_partition(int partition, temp_file_map_t& temp_file_map, const seqset& the_seqset) const
{
	auto temp_file_ptr = std::make_shared<scoped_temp_file>();
	SPLOG("Partition %d: Flattening seqset to temp file \"%s\"."
		, partition, temp_file_ptr->path().c_str());

	dna_sequence start_sequence = find_partition_sequence(partition, m_num_threads);
	seqset_range start_context = the_seqset.find(start_sequence);
	uint64_t start = start_context.begin();
	uint64_t end = the_seqset.size();
	if (partition != (m_num_threads - 1)) {
		dna_sequence end_sequence = find_partition_sequence(partition + 1, m_num_threads);
		seqset_range end_context = the_seqset.find(end_sequence);
		end = end_context.begin();
	}

	SPLOG("Partition %d: start = %lu, end = %lu", partition, start, end);

	file_writer temp_file_writer(temp_file_ptr->path());
	SPLOG("Partition %d entry count: %lu", partition, end - start);
	dna_writer seq_writer(temp_file_writer);
	for(uint64_t i = start; i < end; i++) {
		// atomic progress tracker
		m_cur_progress++;
		dna_sequence the_seq {the_seqset.ctx_entry(i).sequence()};
		if (m_max_read_size > -1 && the_seq.size() > size_t(m_max_read_size)) {
			the_seq = the_seq.subseq(0, m_max_read_size);
		}
		seq_writer.write(the_seq);
		if (((i - start) % 10000000) == 0) {
			SPLOG("Partition %d: Wrote %lu sequences to flat file", partition, i - start);
		}
	}
	SPLOG("Partition %d: Done writing.", partition);

	{{
		std::unique_lock<std::mutex> map_lock(m_temp_file_map_mutex);
		temp_file_map.emplace(partition, temp_file_ptr);
	}}
	SPLOG("Partition %d: Finished flattening.", partition);
}

// Partition DNA space, so map partition and total partition to dna sequence.
// This function takes advantage of the correspondence between base 4 arithmetic
// and DNA sequences.  If the number of threads, and thus the number of partitions
// is a power of 4, then each partition is simply the corresponding integer, i.e.
// a kmer_t.  If it's not a power of 4 (but still a power of 2 as asserted), we
// need to stretch the partitions by a factor of two.
dna_sequence flatten_seqset::find_partition_sequence(int partition, int thread_count)
{
	int total_partitions = thread_count;

	kmer_t partition_kmer = partition;
	if (! is_power_of_4(total_partitions)) {
		partition_kmer *= 2;
		total_partitions *= 2;
	}

	int partition_seq_length = 0;
	while (total_partitions >>= 2) {
		partition_seq_length++;
	}

	return dna_sequence(partition_kmer, partition_seq_length);
}

bool flatten_seqset::is_power_of_2(int n)
{
	if (n <= 0) return false;
	return ! (n & (n - 1));
}

bool flatten_seqset::is_power_of_4(int n)
{
	if (! is_power_of_2(n)) return false;
	return n & 0x55555555;
}

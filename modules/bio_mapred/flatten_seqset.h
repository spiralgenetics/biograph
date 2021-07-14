#pragma once

#include <vector>
#include <string>
#include <map>
#include <memory>
#include <future>

#include "modules/bio_base/dna_sequence.h"

class seqset;
class scoped_temp_file;

// Writes a seqset to a set of temp files that can later be merged into a new seqset.
// Each file covers a power of two divisor of DNA space with a minimum of four files
// (one for each base), thus partitioning the DNA space.  Each partition is assigned
// a thread.
//
// Partitions   Range Covered (in C++ style intervals, i.e. end is not in the interval)
//      4       A-C, C-G, G-T, T-<end>
//      8       AA-AG, AG-CA, CA-CG, CG-GA, GA-GG, GG-TA, TA-TG, TG-<end>
//     16       AA-AC, AC-AG, AG-AT, AT-CA, ... , TA-TC, TC-TG, TG-TT, TT-<end>
//     32       AAA-AAG, AAG-ACA, ACA-ACG, ACG-AGA, ... , TTA-TTG, TTG-<end>
//     64       AAA-AAC, AAC-AAG, AAG-AAT, AAT-ACA, ... , TTA-TTC, TTC-TTG, TTG-TTT, TTT-<end>
class flatten_seqset
{
	using temp_file_map_t = std::multimap<int, std::shared_ptr<scoped_temp_file>>;

public:
	flatten_seqset(std::vector<std::string> seqset_files, int num_threads, int max_read_size = -1);

	std::multimap<int, std::shared_ptr<scoped_temp_file>> operator() ();

	// Return the DNA sequence for the beginning of a given partition in DNA space.
	static dna_sequence find_partition_sequence(int partition, int thread_count);
	static bool is_power_of_2(int n);
	static bool is_power_of_4(int n);

private:
	std::vector<std::string> m_seqset_files;
	int m_num_threads;
	int m_max_read_size;

	mutable std::mutex m_temp_file_map_mutex;
	mutable std::atomic<size_t> m_cur_progress;

	void flatten_partition(int partition, temp_file_map_t& temp_file_map, const seqset& the_seqset) const;
};

// Takes the merged flattened files in DNA space order and produces a seqset
class unflatten_seqset
{
	unflatten_seqset(std::string seqset_file_name, std::vector<std::shared_ptr<scoped_temp_file>>);
};

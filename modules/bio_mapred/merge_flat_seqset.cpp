#include <cassert>
#include <iostream>
#include <memory>
#include <thread>

#include "base/base.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_format/dna_io.h"
#include "modules/bio_mapred/flatten_seqset.h"
#include "modules/bio_mapred/merge_flat_seqset.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/utils.h"
#include "modules/mapred/temp_file.h"

bool check_and_advance(const dna_sequence& main_seq, multi_file_dna_buffer& candidate_buffer);

void merge_flat_seqsets::operator() (
	const std::string& merged_seqset_path
	, const std::multimap<int, std::shared_ptr<scoped_temp_file>>& temp_files_map
	, const bool write_flat, unsigned max_read_len
)
{
	print_progress(0.0);
	std::vector<std::shared_ptr<scoped_temp_file>> merged_temp_files {merge_seqs(temp_files_map)};
	CHECK(flatten_seqset::is_power_of_2(merged_temp_files.size()));
	CHECK(merged_temp_files.size() >= 4);
	create_seqset(merged_seqset_path, merged_temp_files, write_flat, max_read_len);
}

void merge_flat_seqsets::create_seqset(
	const std::string& seqset_path
	, std::vector<std::shared_ptr<scoped_temp_file>> merged_temp_files
	, const bool write_flat,
    unsigned max_read_len
) const
{
	uint64_t sequence_count = m_seq_count.load();
	SPLOG("Creating file builder at \"%s\" with %lu entries.", seqset_path.c_str(), sequence_count);

	std::unique_ptr<spiral_file_create_mmap> builder;
	std::unique_ptr<seqset> build_seqset;
	seqset* the_seqset;
    builder.reset(new spiral_file_create_mmap(seqset_path));
    build_seqset.reset(new seqset(builder->create(), sequence_count, max_read_len));
    the_seqset = build_seqset.get();
	SPLOG("Merging seqset \"%s\" with %lu entries", seqset_path.c_str(), sequence_count);

	// The temp files range over the entire DNA space, so the entire vector of
	// temp files covers all of the future seqset entries.
	multi_file_dna_buffer overall_cursor(merged_temp_files);

	// Each quartile of the temp file vector covers sequences that begin with the
	// respective bases, so we define a cursor into the four subspaces begining with
	// the four DNA bases.
	dna_base_array<boost::optional<multi_file_dna_buffer>> base_cursor;
	CHECK_EQ(0, merged_temp_files.size() % dna_base::k_num_bases)
	    << merged_temp_files.size();
	for (dna_base b : dna_bases()) {
	  base_cursor[b].emplace(temp_files_t(
	      merged_temp_files.begin() +
	      int(b) * merged_temp_files.size() / dna_base::k_num_bases,
	      merged_temp_files.begin() +
	      (int(b) + 1) * merged_temp_files.size() / dna_base::k_num_bases));
	}

	//Flat File 
	std::string flat_out_path = seqset_path + ".flat";
	file_writer flat_out(flat_out_path);
	dna_writer flat_out_dna(flat_out);
	if (! write_flat) {
		flat_out_dna.close();
		flat_out.close();
		std::remove(flat_out_path.c_str());
	} else {
		SPLOG("Writing flat file %s", flat_out_path.c_str());
	}
	
	long entry = 0;
	dna_sequence previous_sequence;
	while (! overall_cursor.at_eof()) {
		const dna_sequence& current_sequence = overall_cursor.get_sequence();
		if (write_flat) {
			//dna_sequence out_seq(cur.begin(), cur.end());
			flat_out_dna.write(current_sequence);
		}
		
		for (dna_base b : dna_bases()) {
		  the_seqset->set_bit(entry, b, check_and_advance(current_sequence, *base_cursor[b]));
		}
                the_seqset->set_entry_size(entry,
                                          uint8_t(current_sequence.size()));
                uint8_t shared = 0;
		for (size_t i = 0; i < std::min(current_sequence.size(), previous_sequence.size()); i++) {
			if (current_sequence[i] != previous_sequence[i]) break;
			shared++;
		}
		the_seqset->set_shared(entry, shared);

		previous_sequence = std::move(current_sequence);
		overall_cursor.advance();
		entry++;

		if ((entry % 400000) == 0) {
			print_progress(float(entry) / float(sequence_count));
		}

		if ((entry % 100000000) == 0) {
			SPLOG("Merged %ld entries of %lu", entry, sequence_count);
		}
	}

	SPLOG("Finalizing seqset merge");
	the_seqset->finalize();
	print_progress(1.0);
	std::cout << std::endl;
}

std::vector<std::shared_ptr<scoped_temp_file>>
	merge_flat_seqsets::merge_seqs
		(const std::multimap<int, std::shared_ptr<scoped_temp_file>>& temp_files_map)
{
	CHECK(! temp_files_map.empty());
	std::vector<std::shared_ptr<scoped_temp_file>> merged_files;

	int partition_count = 0;
	for (auto it = temp_files_map.begin(), end = temp_files_map.end();
			it != end;
			it = temp_files_map.upper_bound(it->first))
	{
		partition_count++;
	}
	merged_files.resize(partition_count);

	int thread_count = 0;
	std::vector<std::future<void>> are_threads_done;
	// Walk unique keys in the multimap.
	for (auto it = temp_files_map.begin(), end = temp_files_map.end();
			it != end;
			it = temp_files_map.upper_bound(it->first))
	{
		int partition = it->first;
		auto thread_lambda = [this, partition, &temp_files_map, &merged_files]() {
			try {
				merged_files[partition] = this->do_merge(partition, temp_files_map);
			} catch (const io_exception& e) {
				SPLOG("Partition %d: Exception \"%s\"", partition, e.message().c_str());
				throw;
			}
		};
		are_threads_done.emplace_back(std::async(std::launch::async, thread_lambda));
		thread_count++;
	}

	for (int i = 0; i < thread_count; i++) {
		are_threads_done[i].get();
	}
	are_threads_done.clear();

	return merged_files;
}


std::shared_ptr<scoped_temp_file> merge_flat_seqsets::do_merge(
	int partition,
	const std::multimap<int, std::shared_ptr<scoped_temp_file>>& temp_files_map
)
{
	using map_value_t = std::multimap<int, std::shared_ptr<scoped_temp_file>>::value_type;

	SPLOG("Partiton: %d Beginning file merge", partition);

	// Create dna_buffers for every value for our partition.
	auto partition_temp_files = temp_files_map.equal_range(partition);
	std::vector<std::unique_ptr<dna_buffer>> dna_buffers;
	std::transform(
		partition_temp_files.first
		, partition_temp_files.second
		, std::back_inserter(dna_buffers)
		, [](map_value_t map_value) { return make_unique<dna_buffer>(map_value.second->path()); }
	);

	// Create the output file and dna_writer.
	auto merge_target_file = std::make_shared<scoped_temp_file>();
	file_writer merge_writer(merge_target_file->path());
	dna_writer merge_target_writer(merge_writer);
	SPLOG("Partition %d: merging to temp file \"%s\"", partition, merge_target_file->path().c_str());

	// Walk the source files looking for the lexicagraphically largest superset of the
	// smallest value.  I.e. AGGA beats AGG, but not AAG
	dna_sequence current_sequence;
	int subseq_count = 0;
	while (!dna_buffers.empty()) {
		current_sequence.clear();
		for (auto dna_buffer_iter = dna_buffers.cbegin()
			; dna_buffer_iter != dna_buffers.cend()
			; dna_buffer_iter++
		) {
			CHECK(*dna_buffer_iter);
			const dna_sequence& it_seq = (*dna_buffer_iter)->get_sequence();
			size_t min_size = std::min(current_sequence.size(), it_seq.size());

			if (current_sequence.subseq(0, min_size) == it_seq.subseq(0, min_size)) {
				if (it_seq.size() > min_size) {
					current_sequence = it_seq;
					subseq_count++;
				}
				continue;
			}
			if (it_seq < current_sequence) {
				current_sequence = it_seq;
			}
		}

		auto dna_buffer_iter = dna_buffers.begin();
		while (dna_buffer_iter != dna_buffers.end()) {
			CHECK(*dna_buffer_iter);
			const dna_sequence& it_seq = (*dna_buffer_iter)->get_sequence();
			while (current_sequence.size() >= it_seq.size()
				&& current_sequence.subseq(0, it_seq.size()) == it_seq) {
				(*dna_buffer_iter)->advance();
				if ((*dna_buffer_iter)->at_eof()) {
					dna_buffer_iter = dna_buffers.erase(dna_buffer_iter);
					goto skip_increment;
				}
			}
			dna_buffer_iter++;
			skip_increment: continue;
		}

		merge_target_writer.write(current_sequence);
		++m_seq_count;
	}

	SPLOG("Partition %d: Finished %lu sequences", partition, m_seq_count.load());
	merge_writer.close();
	return merge_target_file;
}

// Compare the main DNA sequence with an overlap candidate.  To qualify, the
// candidate suffix (i.e. everything but the first base) must be a prefix of
// the main sequence.  If there's an overlap, advance the candidate sequence
// buffer and return true.
bool check_and_advance(const dna_sequence& main_seq, multi_file_dna_buffer& candidate_buffer)
{
	if (candidate_buffer.at_eof()) return false;

	// The potential overlap is the smaller of the candidate suffix or main sequence
	size_t overlap_size = std::min(candidate_buffer.get_sequence().size() - 1, main_seq.size());
//	if (candidate_buffer.get_sequence().subseq(1, overlap_size) == main_seq.subseq(0, overlap_size)) {
	if (subseq_compare(
			candidate_buffer.get_sequence().begin() + 1
			, main_seq.begin()
			, overlap_size
			, overlap_size
                       ) == dna_compare_result::EQUAL) {
		candidate_buffer.advance();
		return true;
	}
	return false;
}


#include <algorithm>
#include <numeric>
#include <atomic>
#include <thread>

#include "base/base.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_format/dna_io.h"
#include "modules/bio_mapred/mem_seqset.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_vector.h"
#include "modules/io/spiral_file.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/io/uuid.h"
#include "modules/io/parallel.h"
#include "modules/mapred/manifest_parallel.h"
#include "modules/mapred/resource_manager.h"
#include "modules/mapred/temp_file.h"

const uint64_t BASES_PER_BYTE = 4;
const uint64_t PADDING_HACK = 32; // goes back to 1 when dna_sequence is refactored
const uint32_t k_max_read_len = 255;  // Maximum read size

REGISTER_TASK(mem_seqset_task);

void mem_seqset_task::validate()
{
	SPLOG_P(LOG_DEBUG,"mem_seqset_task::validate> num_threads: %lu max_mem: %lu, is_paired: %d", num_threads, max_mem, is_paired);

	if (!num_threads) {
		num_threads = 4;
		SPLOG_P(LOG_DEBUG,"mem_seqset_task::validate> threads unspecified, setting to %lu", num_threads);
	}
}

// Given an mmap_vector of a known capacity, make an actual file on disk to hold the data
// Called to before the use of a new mmap_vector
template<class T>
static void reify(mmap_vector<T>& vec) {
	resource_manager rm;
	rm.create_resource(vec.get_buffer(), vec.buffer_size());
	vec.resize(vec.capacity());
}

// Close and truncate a mmap_vector, and return the name of the file which holds its contents
// Called after the use of an mmap_vector
template<class T>
static std::string make_file(mmap_vector<T>& vec) {
	size_t cur_size = vec.size();
	SPLOG("Making file out of vector, size = %lu", cur_size);
	SPLOG("New output size is: %lu, truncating to %lu", cur_size , cur_size * sizeof(T));
	std::string path = vec.get_buffer().path();
	vec.get_buffer().sync();
	vec.get_buffer().close();
	if (::truncate(path.c_str(), cur_size * sizeof(T)) != 0) {
		throw io_exception(std::string("Unable to truncate file: ") + strerror(errno));
	}
	return path;
}

// Convert a flyweight into a dna_slice
static dna_slice decode_flyweight(const char* repo, const flyweight& f) {
	return dna_slice(dna_const_iterator((const unsigned char*) repo, f.start, f.flipped), f.length);
}

// Compare two flyweights as sequences, return true if a is lexicographically before b
static bool flyweight_lt(const char* repo, const flyweight& a, const flyweight& b) {
	if (a.empty) return false;
	if (b.empty) return true;
	dna_compare_result cmp = subseq_compare(
		dna_const_iterator((const unsigned char*) repo, a.start, a.flipped),
		dna_const_iterator((const unsigned char*) repo, b.start, b.flipped),
		a.length, b.length);
	if (cmp != dna_compare_result::EQUAL) {
      return cmp == dna_compare_result::FIRST_IS_LESS ||
          cmp == dna_compare_result::FIRST_IS_PREFIX;
	}
	return *reinterpret_cast<const uint64_t*>(&a) < *reinterpret_cast<const uint64_t*>(&b);
}

// Compare two flyweights as sequences, return true if a is a prefix of b or b is a prefix of a
static bool flyweight_pe(const char* repo, const flyweight& a, const flyweight& b) {
	return subseq_equal(
		dna_const_iterator((const unsigned char*) repo, a.start, a.flipped),
		dna_const_iterator((const unsigned char*) repo, b.start, b.flipped),
		std::min(a.length, b.length));
}

// Very similar but not quite std::unique
// Basically, std::unique makes the output the 'first' element of the group, I need
// the 'last' element of a group (since I'm using pe as my predicate, I want the longest
// prefix).  Maybe I could fix this with reverse iterators or something, but whatever
template<class It, class Func>
It my_unique(const It& begin, const It& end, const Func& f) {
	It cur_out = begin;
	It cur_in = begin;
	while(cur_in != end) {
		*cur_out = *cur_in;
		cur_in++;
		if (cur_in == end) { cur_out++; break; }
		if (!f(*cur_out, *cur_in)) {
			cur_out++;
		}
	}
	return cur_out;
}

// Expand one read: i.e., compute the set of all suffixes, until we find a suffix that is in the
// list of originals. For the read, ACTGGA for example, we would add CTGGA, TGGA, etc, until perhaps
// GGA was the prefix of another sequence in 'originals', say GGATTA. We put the resulting outputs
// as flyweights into output (presumed empty initially).
// We require the repo (to decode flyweights), and the list of originals.
size_t mem_seqset_task::expand_one_read(tracked_vector<flyweight>& output, flyweight read)
{
	// Compute size and direction
	size_t read_len = read.length;
	int sign = 1 - 2*read.flipped;
	// Go over each offset
	for(size_t offset = 1; offset < read_len; offset++) {
		// Compute new flyweight
		size_t sz = read_len - offset;
		flyweight f2(read.start + sign * offset, sz, read.flipped);

		// See if the new sequence is a prefix of an original
		bool is_prefix = std::binary_search(m_originals->begin(), m_originals->end(), f2,
			[this, sz](const flyweight& a, const flyweight& b) {
				return decode_flyweight(m_repo, a).subseq(0, std::min(sz, size_t(a.length))) <
				decode_flyweight(m_repo, b).subseq(0, std::min(sz, size_t(b.length)));
		});
		// If so, break
		if (is_prefix) {
			return read_len - offset;
		}

		// Push back to output
		output.push_back(f2);
	}
	return 0;
}

// One pass performs an expansion 'originals' into a new output memory map constructed locally and
// and returned as a file for later merge sorting.  This expansion is multithreaded, using atomics.
// During this process, we move 'next_read' forward, up until we hit the end of originals, which we may
// or may not do in a single pass.  We also keep track of the highest lexicographical subsequence we ever
// generate (as 'worst_ever')
std::string mem_seqset_task::one_expand_pass(const progress_handler_t& progress)
{
	CHECK_GT(m_originals->size(), 0);
	// Since each thread only checks for stopping once per read, we need to stop early, lest each thread
	// sees that there is space for it, but there is not enough space for *all* threads at once
	uint64_t high_water = m_max_buf_size - num_threads * k_max_read_len;
	CHECK_GT(high_water, 0);
	// Construct output mmap
	mmap_vector<flyweight> out(m_max_buf_size);
	reify(out);
	// Make an atomic counter to maintain output position
	std::atomic<uint64_t> next_write(0);
	std::vector<std::thread> threads;
	SPLOG_P(LOG_DEBUG, "Starting threads");
	// Kick off all our threads
	for(size_t j = 0; j < num_threads; j++) {
		threads.emplace_back([&]() {
			// Init estimated output position to '0'
			uint64_t maxo = 0;
			// Vector for outputs
			tracked_vector<flyweight> output(track_alloc("mem_seqset:flyweight"));
			// While our next output position is below high-water
			while(maxo < high_water) {
				// Get next read to process, since multiple threads are pulling
				// we may get an entry past the end, in which case we leave the loop
				uint64_t i = m_next_read.fetch_add(1, std::memory_order_relaxed);
				if (i >= m_originals->size()) break;
				if (i % (1*1024*1024) == 0) {
					SPLOG_P(LOG_DEBUG, "Processing entry %lu", i);
				}
				// Grab the data and process
				expand_one_read(output, (*m_originals)[i]);
				// If nothing to report, forget it
				if (output.size() == 0) continue;
				// Atomically allocate some space in the shared ouput area
				uint64_t o = next_write.fetch_add(output.size(), std::memory_order_relaxed);
				// Copy the data over
				for(size_t k=0; k < output.size(); k++) {
					out[o++] = output[k];
				}
				// Clear output vector for next go around
				output.clear();
				maxo = o;
			}
		});
	}
	while(true) {
		uint64_t pos = m_next_read.load();
		uint64_t opos = next_write.load();
		if (pos >= m_originals->size() || opos >= high_water) break;
		progress(double(pos) / double(m_originals->size()));
		sleep(1);
	}
	// Await all the results
	for(size_t i = 0; i < num_threads; i++) {
		threads[i].join();
	}
	uint64_t out_size = next_write.load();
	CHECK_GT(out_size, 0);
	double cur_progress = double(m_next_read.load()) / double(m_originals->size());
	std::string out_file;

	lambda_watchdog(progress, cur_progress, [&]() {
		// Get the final size of the output, should be somewhere between high_water and k_hard_max
		// Resize output to match
		SPLOG("Output size is: %lu", out_size);
		out.resize(out_size);
		// Sort all of the ouputs
		SPLOG("Doing out sort");
		parallel_sort_in_place(out.begin(), out.end(),
			[this](const flyweight& a, const flyweight& b) { return flyweight_lt(m_repo, a, b); });
		// Unique the outputs

		SPLOG("Dedupping out");
		auto new_end2 = my_unique(out.begin(), out.end(),
			[this](const flyweight& a, const flyweight& b) { return flyweight_pe(m_repo, a, b); });
		out.resize(new_end2 - out.begin());

		SPLOG("Output size is: %lu", out.size());
		// Update worst_ever
		if (decode_flyweight(m_repo, out[out.size() - 1]) > decode_flyweight(m_repo, m_worst_ever)) {
			m_worst_ever = out[out.size() - 1];
		}
		// Finalize the physical output file
		out_file = make_file(out);
	});
	return out_file;
}

// Merge multiple files into a single file, and compute the location of the first entry with A, C,
// G, and T (a.k.a the C(a) table). Writes the results into a single large file.  Operates in
// multiple passes, since we use *sort* to do the merge since it's faster to use parallel sort than
// to do serial finger-merge amazingly enough. We should probably use a parallel merge, but it was
// complex, and not the bottleneck. Basically, we presume each file is evenly distributed, load a
// chunk of each into a buffer, sort the buffer and then copy the results out, moving each file
// forward by however much was safely consumed (since our estimates are not certain).
size_t mem_seqset_task::do_merge(file_writer& output, std::vector<file_reader>& inputs)
{
	m_base_pos[0] = 0;
	// Put remaining size of each file into a vector, and compute the total remaining
	std::vector<size_t> remaining(inputs.size());
	size_t tot_remaining = 0;
	for(size_t i = 0; i < inputs.size(); i++) {
		remaining[i] = inputs[i].size() / sizeof(flyweight);
		tot_remaining += remaining[i];
	}
	// Make my merge buffer
	tracked_vector<flyweight> merge_buf(m_max_buf_size, track_alloc("mem_seqset::do_merge:flyweight"));
	// While there is anything remaining
	size_t output_size = 0;
	while(tot_remaining) {
		// If this is not the first run, we have a single entry from the last run
		// This entry is there to make sure prefix de-dup operates correctly
		size_t buf_off = 0;
		if (output_size != 0) {
			buf_off = 1; // Hold over starts in the bucket
		}
		// Compute where to load each group of sequences by estimating the
		// proportional size of each file on the assumption of even distribution
		// which should be close to true
		std::vector<size_t> buf_start(inputs.size());
		std::vector<size_t> estimates(inputs.size());
		for(size_t i = 0; i < inputs.size(); i++) {
			// New buffer starts where the last one left off
			buf_start[i] = buf_off;
			// All of the entries will fit in out buffer, don't bother estimating
			// Note: tot_remaining must be strictly less, since we may have one hold-over entry
			if (tot_remaining < m_max_buf_size) {
				estimates[i] = remaining[i];
			} else {
				// Otherwise, presume proportional useage
				estimates[i] = remaining[i] * (m_max_buf_size- inputs.size() - 1) / tot_remaining + 1;
			}
			// Add used amount to current buffer offset
			buf_off += estimates[i];
		}
		SPLOG("Total size = %lu", buf_off);
		// Read all the data into the buffer.  While we are doing that we also
		// check each inputs 'highest' value, and find the 'lowest' of any input
		// Any entries larger than this must be ignored, since the valid area of the merge
		// is the region over which all files have entries.  We begin be setting our 'min' to
		// the worst case 'max', which we found during our earlier run
		flyweight lowest = m_worst_ever;
		for(size_t i = 0; i < inputs.size(); i++) {
			SPLOG("i = %d, reading %lu entries", (int) i, estimates[i]);
			// If I put no entries from this input into the mix, ignore it
			if (estimates[i] == 0) continue;
			// Physically load the data based on prior calculations
			size_t r = inputs[i].read((char *) &merge_buf[buf_start[i]], estimates[i] * sizeof(flyweight));
			if (r != estimates[i] * sizeof(flyweight)) {
				throw io_exception("Incomplete read");
			}
			// If we will fit all data, we are valid, no need to compute lowest
			if (tot_remaining < m_max_buf_size) {
				continue;
			}
			// Find the last entry for this file
			flyweight end_of_file = merge_buf[buf_start[i] + estimates[i] - 1];
			// Update 'lowest'
			if (decode_flyweight(m_repo, end_of_file) < decode_flyweight(m_repo, lowest)) {
				lowest = end_of_file;
			}
			SPLOG("last = %s", decode_flyweight(m_repo, end_of_file).as_string().c_str());
		}
		SPLOG("Lowest EOF = %s", decode_flyweight(m_repo, lowest).as_string().c_str());
		// Find the last valid entry in each buffer and compact
		// IE: remove any entries > lowest
		size_t final_size = 0;
		// Keep hold-over regardless
		if (output_size != 0) {
			final_size = 1; // Hold over starts in the bucket
		}
		// For each input
		for(size_t i = 0; i < inputs.size(); i++) {
			// Compute useful size of this input by searching for 'lowest' in sorted region
			size_t actual_size = std::upper_bound(
				merge_buf.begin() + buf_start[i],
				merge_buf.begin() + buf_start[i] + estimates[i],
				lowest,
				[this](const flyweight& a, const flyweight& b) {
					return decode_flyweight(m_repo, a) < decode_flyweight(m_repo, b);
				}) - (merge_buf.begin() + buf_start[i]);
			SPLOG("i = %d, actual size = %lu", (int) i, actual_size);
			// We will consume this much, update remaining for next time around
			remaining[i] -= actual_size;
			tot_remaining -= actual_size;
			// Adjust our file point to be ready to read from here
			inputs[i].seek(inputs[i].size() - (remaining[i] * sizeof(flyweight)));
			// Move the data down, compacting away data that will be ignored, and track total valid size
			memmove(&merge_buf[final_size], &merge_buf[buf_start[i]], actual_size * sizeof(flyweight));
			final_size += actual_size;
		}
		SPLOG("Sorting top %lu entries", final_size);
		// Do the actual parallel sort all of this set the stage for
		parallel_sort_in_place(merge_buf.begin(), merge_buf.begin() + final_size,
			[this](const flyweight& a, const flyweight& b) { return flyweight_lt(m_repo, a, b); });
		// Unique the sorted data
		SPLOG("Uniqing top %lu entries", final_size);
		auto deduped_it = my_unique(merge_buf.begin(), merge_buf.begin() + final_size,
			[this](const flyweight& a, const flyweight& b) { return flyweight_pe(m_repo, a, b); });
		size_t uniq_size = deduped_it - merge_buf.begin();
		SPLOG("Final output size = %lu", uniq_size);
		// Now we compute the location of a top level break (IE: ATTTTC - CAAAGG or whatever), basically
		// where in the final sorted order we first see any entry with a new initial base pair
		// To do this, we check the last value matches the first, but of course all 3 breaks may be in the
		// same buffer, so we need to account for that
		size_t letter_x = uniq_size - 1;
		while (decode_flyweight(m_repo, merge_buf[0])[0] != decode_flyweight(m_repo, merge_buf[letter_x])[0]) {
			// Have a break, it's for 'base'
			int base = (int) decode_flyweight(m_repo, merge_buf[letter_x])[0];
			SPLOG("Computing letter break for base %d (%c)", base, (char) dna_base(base));
			// Make a single base-pair flyweight
			flyweight letter = merge_buf[letter_x];
			letter.length = 1;
			// Seek for its lower bound in the region under consideration
			auto it = std::lower_bound(merge_buf.begin(), merge_buf.begin() + letter_x + 1, letter,
				[this](const flyweight& a, const flyweight& b) { return flyweight_lt(m_repo, a, b); });
			// Store the results
			size_t offset = it - merge_buf.begin();
			m_base_pos[base] = output_size + offset;
			SPLOG("Result is %lu, seq = %s", m_base_pos[base],
				decode_flyweight(m_repo, merge_buf[offset]).as_string().c_str());
			// Move before that point to see if there are other breaks to process
			letter_x = offset - 1;
			SPLOG("Prev seq = %s",
				decode_flyweight(m_repo, merge_buf[letter_x]).as_string().c_str());
		}
		// If we have any more to do
		if (tot_remaining) {
			// Write everything but the last entry, which we hold over in position 0
			output.write((const char*) &merge_buf[0], (uniq_size - 1) * sizeof(flyweight));
			output_size += uniq_size - 1;
			merge_buf[0] = merge_buf[uniq_size - 1];
		} else {
			// Write it all
			output.write((const char*) &merge_buf[0], uniq_size * sizeof(flyweight));
			output_size += uniq_size;
		}
	}
	// Close the merged output
	output.close();
	// Clear any inputs, allowing files to close
	inputs.clear();
	// Put a sentinal into the base_pos array
	m_base_pos[4] = output_size;
	// Return the output size
	return output_size;
}

// fquery represents a logical 'cursor' into a file of flyweights from 'start' till 'end'
// which keeps a 1M buffer for the physical data.  Basically we need 5 of these (one for
// each base pair, and one the whole file) to construct the seqset
class fquery
{
public:
	fquery(const char* repo, const std::string& fly_file, size_t start, size_t end)
		: m_repo(repo)
		, m_file(fly_file)
		, m_buf(1024*1024, track_alloc("mem_seqset:fquery"))
		, m_buf_pos(0)
		, m_buf_end(0)
		, m_cur(start)
		, m_end(end)
	{
		// Go to start and read initial buffer
		m_file.seek(sizeof(flyweight) * start);
		read_next();
	}
	bool has_more() { return m_cur != m_end; }
	size_t get_buf_pos() const { return m_buf_pos; }
	flyweight get_flywt() const { return m_buf[m_buf_pos]; }
	dna_slice current() { return decode_flyweight(m_repo, m_buf[m_buf_pos]); }
	void next() {
		// Move cur and buffer forward, refill buffer as needed
		m_cur++;
		m_buf_pos++;
		if (m_buf_pos == m_buf_end) {
			read_next();
		}
	}
	// Compare the current position to see the main entry is a near suffix of the current
	// position, that is, for main of CTAGG, TCTAG would qualify, but TCTAA would not
	bool check_move(const dna_slice& main) {
		// If no more entries, forget it
		if (!has_more()) return false;
		// Check for near suffix condition
		flyweight f = m_buf[m_buf_pos];
		int size = std::min(f.length - 1, (int) main.size());
		if (decode_flyweight(m_repo, f).subseq(1, size) == main.subseq(0, size)) {
			// If so, we also more forward one
			next();
			return true;
		}
		// We should never get 'behind' the main entry, blow up if we do
		if (decode_flyweight(m_repo, f).subseq(1, f.length - 1) < main) {
			SPLOG("Inconsistency in seqset generation!");
			SPLOG("%s", decode_flyweight(m_repo, f).as_string().c_str());
			SPLOG(" %s", main.as_string().c_str());
			throw io_exception("Inconsistency in seqset generation!");
		}
		return false;
	}
private:
	void read_next() {
		// If no more to read, stop
		if (m_cur == m_end) return;
		// Pick buffer size + read
		size_t size = std::min(m_buf.size(), m_end - m_cur);
		size_t r = m_file.read((char *) &m_buf[0], sizeof(flyweight) * size);
		if (r != sizeof(flyweight) * size) {
			throw io_exception("Incomplete fquery read");
		}
		// Update buffer data
		m_buf_pos = 0;
		m_buf_end = size;
	}
	const char* m_repo;  // The repo (to decode flyweights)
	file_reader m_file;  // THe file we are reading from
	tracked_vector<flyweight> m_buf;  // The current buffer
	size_t m_buf_pos;  // Our position in the buffer
	size_t m_buf_end;  // The end fo the buffer
	size_t m_cur;  // Our logical current position
	size_t m_end;  // Our logical end
};

manifest mem_seqset_task::output_seqset(const progress_handler_t& progress, const std::string& output_name, size_t tot_size, bool write_flat)
{
	SPLOG("Running final step");
	std::string mmap_file{ ::path(CONF_S(path_bulkdata)).bare_path() + "/" + make_uuid() };

	std::unique_ptr<spiral_file_create_mmap> builder;
	std::unique_ptr<seqset> new_seqset;
	seqset* the_seqset;

    builder.reset(new spiral_file_create_mmap(mmap_file));
    new_seqset.reset(new seqset(builder->create(), tot_size, k_max_read_len));
    the_seqset = new_seqset.get();

	// Make queries for each base and for main
	SPLOG("tot_size = %lu, m_base_pos[0] = %lu, m_base_pos[1] = %lu, m_base_pos[2] = %lu, m_base_pos[3] = %lu, m_base_pos[4] = %lu"
		, tot_size, m_base_pos[0], m_base_pos[1], m_base_pos[2], m_base_pos[3], m_base_pos[4]);
	fquery qm(m_repo, output_name, 0, tot_size);
	dna_base_array<boost::optional<fquery>> qbase;
	for (dna_base b : dna_bases()) {
	  qbase[b].emplace(m_repo, output_name, m_base_pos[int(b)], m_base_pos[int(b)+1]);
	}
	dna_slice prev;
	// Go over each entry
	subprogress sp_output(progress, 0.00, 0.85);

	//Flat File
	std::string flat_out_path = mmap_file + ".flat";
	file_writer flat_out(flat_out_path);
	dna_writer flat_out_dna(flat_out);
	if (! write_flat) {
		flat_out_dna.close();
		flat_out.close();
		std::remove(flat_out_path.c_str());
	} else {
		SPLOG("Writing flat file %s", flat_out_path.c_str());
	}

	size_t row = 0;
	while(qm.has_more()) {
		if (row % 1000000 == 0) {
			SPLOG_P(LOG_DEBUG, "Doing row %lu", row);
		}
		dna_slice cur = qm.current();

		if (write_flat) {
			dna_sequence out_seq(cur.begin(), cur.end());
			flat_out_dna.write(out_seq);
		}

		// Set the appropriate bits
        for (dna_base b: dna_bases()) {
          the_seqset->set_bit(row, b, qbase[b]->check_move(cur));
        }
                the_seqset->set_entry_size(row, cur.size());
                unsigned shared = 0;
		// Compute shared overlap
		// TODO: This is SLOW, but it's not in the critical path of the whole process
		for (size_t i = 0; i < std::min(cur.size(), prev.size()); i++) {
				if (cur[i] != prev[i]) break;
				shared++;
		}
		the_seqset->set_shared(row, shared);
		row++;
		sp_output(double(row) / double(tot_size));
		qm.next();
		prev = cur;
	}
	// Finish things up
	SPLOG("Finalizing");

	the_seqset->finalize(subprogress(progress, 0.85, 0.99));

	SPLOG("Writing manifest!");
	// Write as a resource
	manifest mout;

	file_info fi;
	fi.file = mmap_file;
    fi.size = builder->close();
	fi.num_records = 0;
	mout.add(fi, 0);

	return mout;
}

// Main function to generate a new seqset using in memory parallel sorts
manifest mem_seqset_task::do_mem_seqset(const progress_handler_t& progress)
{
	SPLOG("Size of flyweight = %d", (int) sizeof(flyweight));
	// Sort + dedup
	SPLOG("Doing original sort");
	lambda_watchdog(progress, 0.01, [&]() {
		parallel_sort_in_place(m_originals->begin(), m_originals->end(),
			[this](const flyweight& a, const flyweight& b) { return flyweight_lt(m_repo, a, b); });
	});

	SPLOG("Deduping");
	auto last_empty_flywt = m_originals->end();
	lambda_watchdog(progress, 0.03, [&]() {
		while (--last_empty_flywt > m_originals->begin() && last_empty_flywt->empty) {}
		last_empty_flywt++;
	});
	flyweight* new_end = nullptr;
	lambda_watchdog(progress, 0.08, [&]() {
		new_end = my_unique(m_originals->begin(), last_empty_flywt,
			[this](const flyweight& a, const flyweight& b) { return flyweight_pe(m_repo, a, b); });
	});
	m_originals->resize(new_end - m_originals->begin());

	// Set inital 'worst_ever'
	m_worst_ever = (*m_originals)[m_originals->size() - 1];
	SPLOG("New size = %d", (int) m_originals->size());

	if (m_originals->size() * k_max_read_len < m_max_buf_size) {
		m_max_buf_size = (m_originals->size() + num_threads) * k_max_read_len;
	}

	SPLOG("Doing read expansion");
	// Expand all reads, make an atomic to keep track of where we are
	m_next_read.store(0);
	std::vector<file_reader> inputs;
	// While we are not complete
	subprogress sp_expand(progress, 0.15, 0.60);
	while(m_next_read.load() < m_originals->size()) {
		// Do a parallel expansion, and add output to be merged later
		std::string input = one_expand_pass(sp_expand);
		SPLOG("File name = %s", input.c_str());
		inputs.emplace_back(input);
	}
	// Write out 'originals' which also need to be the merge, and add them
	SPLOG("Closing originals mmap");
	std::string orig_file = make_file(*m_originals);
	inputs.emplace_back(orig_file);
	// Run the merge process
	SPLOG("Merging all things");
	scoped_temp_file flat_weight_file;
	std::string output_name{flat_weight_file.path()};
	file_writer output(output_name);
	size_t tot_size = 0;
	lambda_watchdog(progress, 0.65, [&]() {
		tot_size = do_merge(output, inputs);
	});

	return output_seqset(subprogress(progress, 0.70, 1.0), output_name, tot_size, m_write_flat);
}

struct map_ref_reads
{
	map_ref_reads(
		const mem_seqset_task& the_task
		, const reference& ref
		, std::shared_ptr<mmap_vector<flyweight>> the_flyweights
		, int file_info_count
	)
		: m_task(the_task), m_ref(ref), m_base_count(file_info_count), m_flywts_ptr(the_flyweights)
	{
	}

	map_ref_reads(map_ref_reads&& rhs)
		: m_task(rhs.m_task)
		, m_ref(rhs.m_ref)
		, m_base_count(std::move(rhs.m_base_count))
		, m_flywts_ptr(rhs.m_flywts_ptr)
	{
	}

	map_ref_reads& operator=(map_ref_reads&& rhs) = delete;

	map_ref_reads(const map_ref_reads& rhs) = delete;
	map_ref_reads& operator=(const map_ref_reads& rhs) = delete;

	~map_ref_reads()
	{}

	void operator()(
		const std::string& read_id
		, const corrected_reads& the_reads
		, int file_info_id
		, uint64_t record_id
	);

	std::vector<uint64_t> get_base_counts() const { return m_base_count; }

private:
	const mem_seqset_task& m_task;
	const reference& m_ref;
	std::vector<uint64_t> m_base_count; // Per file-info non-reference base count
	std::shared_ptr<mmap_vector<flyweight>> m_flywts_ptr;
};

inline void map_ref_reads::operator()(
	const std::string& read_id
	, const corrected_reads& the_reads
	, int file_info_id
	, uint64_t record_id
)
{
	if (the_reads.empty()) {
		SPLOG("Corrected read ID \"%s\" is completely empty.", read_id.c_str());
		return;
	}

	flyweight the_flyweight;
	dna_sequence the_read = the_reads[0].corrected;
	auto bwt = m_ref.get_bwt();

	auto fwd = bwt.find(the_read);
	// forward match?
	if (fwd.matches() > 0) {
		the_flyweight = flyweight(fwd.get_match(0), the_read.size(), false);
	} else {
		auto rev = bwt.find(the_read.rev_comp());
		// reverse complement match?
		if(rev.matches() > 0) {
			// BWT matches in the forward direction only, but this read's reverse complement matched.
			// Make a flyweight in the forward direction, but track its complement.
			the_flyweight = flyweight(rev.get_match(0), the_read.size(), false).rev_comp();
		} else {
			// Non-reference read
			the_flyweight.flipped = false;
			the_flyweight.non_ref = true;
			m_base_count[file_info_id] += the_read.size();
		}
	}

	(*m_flywts_ptr)[m_task.flywt_index(record_id)] = the_flyweight;
	(*m_flywts_ptr)[m_task.flywt_index(record_id) + 1] = the_flyweight.rev_comp();

	if (the_reads.size() == 2) {
		if (! m_task.is_paired) {
			throw io_exception(printstring(
				"map_ref_reads::operator()> Paired reads found in unpaired task, read id \"%s\"", read_id.c_str()));
		}

		flyweight mate_flyweight;
		dna_sequence mate_read = the_reads[1].corrected;

		auto mate_fwd = bwt.find(mate_read);

		// forward match?
		if (mate_fwd.matches() > 0) {
			mate_flyweight = flyweight(mate_fwd.get_match(0), mate_read.size(), false);
		} else {
			auto mate_rev = bwt.find(mate_read.rev_comp());
			// reverse complement match?
			if (mate_rev.matches() > 0) {
				// BWT matches in the forward direction only, but this read's reverse complement matched.
				// Make a flyweight in the forward direction, but track its complement.
				mate_flyweight = flyweight(mate_rev.get_match(0), mate_read.size(), false).rev_comp();
			} else {
				// Non-reference read
				mate_flyweight.flipped = false;
				mate_flyweight.non_ref = true;
				m_base_count[file_info_id] += mate_read.size();
			}
		}
		(*m_flywts_ptr)[m_task.flywt_index(record_id) + 2] = mate_flyweight;
		(*m_flywts_ptr)[m_task.flywt_index(record_id) + 3] = mate_flyweight.rev_comp();

	} else if (the_reads.size() == 1) {
		if (m_task.is_paired) {
			flyweight empty_flywt;
			empty_flywt.empty = true;
			(*m_flywts_ptr)[m_task.flywt_index(record_id) + 2] = empty_flywt;
			(*m_flywts_ptr)[m_task.flywt_index(record_id) + 3] = empty_flywt;
		}
	} else {
		throw io_exception(printstring("Corrected read ID \"%s\" has %lu mates!  Expected one or two."
			, read_id.c_str(), the_reads.size()));
	}
}

// Used to check the validitiy of all non-empty flyweights. Only used when run_tests == true.
struct validate_flyweights
{
	validate_flyweights(
		const mem_seqset_task& the_task
		, const reference& ref
		, std::shared_ptr<mmap_vector<flyweight>> the_flyweights
		, const char* repo
	)
		: m_task(the_task), m_ref(ref), m_flywts_ptr(the_flyweights), m_repo(repo)
	{}

	void check_one_flyweight(dna_sequence comparison_seq, flyweight fly) const {

		// Every flyweight should be in the repo
		CHECK(fly.length == (decode_flyweight(m_repo, fly).as_string().length()));
		if(fly.non_ref) { return; }

		// Reference flyweights should be in the reference
		dna_const_iterator ref_seq_iter;

		if(fly.flipped) {
			ref_seq_iter = m_ref.get_dna(fly.start).rev_comp();
		} else {
			ref_seq_iter = m_ref.get_dna(fly.start);
		}

		// SPLOG("comparison_seq: %s", comparison_seq.as_string().c_str());
		// SPLOG("  ref_seq_iter: %s", dna_slice(ref_seq_iter, comparison_seq.size()).as_string().c_str());
		CHECK(subseq_equal(ref_seq_iter, comparison_seq.begin(), fly.length));
	}

	void check_one_flyweight(dna_slice slice, flyweight fly) const {
		check_one_flyweight(dna_sequence(slice.as_string()), fly);
	}

	void operator()(
		const std::string& read_id
		, const corrected_reads& the_reads
		, int file_info_id
		, int record_id
	) const
	{
		if (the_reads.size() > 2) {
			throw io_exception(printstring(
				"Unexpected number of reads: read ID \"%s\" has size %lu", read_id.c_str(), the_reads.size()));
		}
		int base_id = m_task.flywt_index(record_id);
		dna_sequence seq = the_reads[0].corrected;

		check_one_flyweight(seq, (*m_flywts_ptr)[base_id]);
		check_one_flyweight(seq.rev_comp(), (*m_flywts_ptr)[base_id + 1]);

		// If paired, check the mate
		if (the_reads.size() == 2) {
			CHECK(m_task.is_paired);
			dna_sequence mate_seq = the_reads[1].corrected;

			check_one_flyweight(mate_seq, (*m_flywts_ptr)[base_id + 2]);
			check_one_flyweight(mate_seq.rev_comp(), (*m_flywts_ptr)[base_id + 3]);
		}
	}

private:
	const mem_seqset_task& m_task;
	const reference& m_ref;
	std::shared_ptr<mmap_vector<flyweight>> m_flywts_ptr;
	const char* m_repo;
};

struct map_non_ref_reads
{
	map_non_ref_reads(
		const mem_seqset_task& the_task
		, std::shared_ptr<mmap_vector<flyweight>> the_flyweights
		, const std::vector<uint64_t>& base_counts
		, char* repo_start
		, size_t ref_size
	)
		: m_task(the_task)
		, m_flywts_ptr(the_flyweights)
		, m_repo_non_ref_offset(ref_size)
	{
		m_repo_offsets.reserve(base_counts.size());
		m_repo_offsets.push_back(0);
		std::partial_sum(base_counts.cbegin(), base_counts.cend(), std::back_inserter(m_repo_offsets));
		m_repo_iter = dna_iterator(reinterpret_cast<unsigned char*>(repo_start), ref_size, false);
		SPLOG("m_repo_iter starts at %p", m_repo_iter.get_data() + m_repo_iter.get_offset());
		SPLOG("m_repo_iter offset = %ld", m_repo_iter.get_offset());
		SPLOG("m_repo_offsets.back() = %lu", m_repo_offsets.back());
	}

	void operator()(
		const std::string& read_id
		, const corrected_reads& the_reads
		, int file_info_id
		, size_t record_id
	);

private:
	const mem_seqset_task& m_task;
	std::shared_ptr<mmap_vector<flyweight>> m_flywts_ptr;
	std::vector<uint64_t> m_repo_offsets; // Each file_info has its own repo area, these are the offsets into them.
	dna_iterator m_repo_iter;
	size_t m_repo_non_ref_offset;

	void fill_non_ref_flywt(flyweight& the_flywt, const dna_sequence& the_read_sequence, int file_info_id);
};

inline void map_non_ref_reads::operator()(
	const std::string& read_id
	, const corrected_reads& the_reads
	, int file_info_id
	, size_t record_id
)
{
	auto flywt_id = m_task.flywt_index(record_id);

	auto& the_flywt = (*m_flywts_ptr)[flywt_id];
	if (the_flywt.non_ref) {
		fill_non_ref_flywt(the_flywt, the_reads[0].corrected, file_info_id);
		(*m_flywts_ptr)[flywt_id + 1] = the_flywt.rev_comp();
	}

	if (m_task.is_paired) {
		auto& mate_flywt = (*m_flywts_ptr)[flywt_id + 2];
		if (mate_flywt.non_ref) {
			fill_non_ref_flywt(mate_flywt, the_reads[1].corrected, file_info_id);
			(*m_flywts_ptr)[flywt_id + 3] = mate_flywt.rev_comp();
		}
	}
}

inline void map_non_ref_reads::fill_non_ref_flywt(flyweight& the_flywt, const dna_sequence& the_read_sequence, int file_info_id)
{
	the_flywt.start = m_repo_non_ref_offset + m_repo_offsets[file_info_id];
	the_flywt.length = the_read_sequence.size();

	// Storing to m_repo_iter in this way is not threadsafe, so hide
	// it behind a lock until DEVEL-370 is fixed.  Note that this
	// effectively removes parallelism for this step, so is not an
	// ideal long term fix.
	static std::mutex mutex;
	std::lock_guard<std::mutex> lock(mutex);

	for (const auto& base : the_read_sequence) {
		*(m_repo_iter + m_repo_offsets[file_info_id]++) = base;
	}
}

void mem_seqset_task::load_repo(const reference& ref, const progress_handler_t& progress)
{
	resource_manager rm;
	manifest_reader mr(input);
	size_t record_count = input.get_num_records();

	SPLOG("Loading reads into in-memory repo");
	SPLOG("Reference size = %lu", ref.size());
	SPLOG("%lu corrected reads bases.", input.metadata().get<size_t>(meta::ns::readonly, "corrected_read_bases"));

	m_originals = std::make_shared<mmap_vector<flyweight>>(count_flyweights());
	reify(*m_originals);

	// Pass 1 through the data.  Make flyweights for reference reads and place-holders
	// for non-ref reads and mates of reads with no mate pair (non-existent reads).
	std::string key;
	corrected_reads value;
	// pre-load so we don't get a race condition later
	auto bwt = ref.get_bwt();
	CHECK(bwt.valid());

	SPLOG("Starting reference read mapping");
	auto functor = map_ref_reads(*this, ref, m_originals, input.count_file_infos());
	std::vector<uint64_t> base_counts;
	lambda_watchdog(progress, 0.4, [&]() {
		base_counts
			= manifest_parallelize<map_ref_reads, std::string, corrected_reads>(input, std::move(functor), subprogress(progress, 0, 0.4)).get_base_counts();
	});
	SPLOG("Reference read mapping complete");
	SPLOG("base_counts size = %lu, base_counts[0] = %lu", base_counts.size(), base_counts[0]);

	size_t accumulated_bases = std::accumulate(base_counts.begin(), base_counts.end(), 0ULL);
	SPLOG("Accumulated bases = %lu", accumulated_bases);
	size_t repo_base_count = ref.size() + accumulated_bases;
	size_t repo_mem_needed = ((repo_base_count + 7) / BASES_PER_BYTE) + PADDING_HACK;
	rm.create_resource(m_repo_mmap, repo_mem_needed);
	SPLOG("Repo memmap buffer starts at %p", m_repo_mmap.buffer());
	std::memcpy(m_repo_mmap.buffer(), ref.get_dna(0).get_data(), bases_to_data_size(ref.size()));
	SPLOG("Record count %zu, base count = %zu, paired = %d", record_count, repo_base_count, is_paired);

	size_t mem_available = max_mem * 1024 * 1024 * 1024;
	// Account for overhead in gnu parallel sort (+40%)
	size_t mem_needed = size_t((repo_mem_needed + m_originals->size() * sizeof(flyweight)) * 1.4);
	if(mem_needed > mem_available) {
		throw io_exception(printstring("Insufficient memory (need %lu bytes, only %lu available)", mem_needed, mem_available));
	}
	mem_available -= mem_needed;
	SPLOG("mem_needed: %lu mem_available: %lu", mem_needed, mem_available);
	m_max_buf_size = size_t(mem_available / sizeof(flyweight));

	// Pass 2 through the data.  Allocate space for the repo now that we know how many reads
	// there are exactly, then build the repo and fill in the non-reference reads.
	SPLOG("Mapping non-ref reads");
	lambda_watchdog(progress, 0.6, [&]() {
		manifest_parallelize<map_non_ref_reads, std::string, corrected_reads>(
			input
			, map_non_ref_reads(
				*this
				, m_originals
				, base_counts
				, m_repo_mmap.buffer()
				, ref.size()
                                ),
            subprogress(progress, 0.4, 0.6)
		);
	});
	SPLOG("Non-reference read mapping complete");
	m_repo = m_repo_mmap.buffer();

	if (run_tests) {
		size_t count = 0;
		size_t empty = 0;
		for (auto flywt : *m_originals) {
			if (flywt.empty) {
				empty++;
				continue;
			}
			if (flywt.non_ref) {
				CHECK(flywt.start >= ref.size());
			} else {
				CHECK(flywt.start < ref.size());
			}
			count++;
		}
		SPLOG("Flyweight preflight test complete! %lu tested, %lu empty flyweights skipped", count, empty);

		SPLOG("Validating flyweights");
		manifest_parallelize<validate_flyweights, std::string, corrected_reads>(input, validate_flyweights(*this, ref, m_originals, m_repo), subprogress(progress, 0.6, 1));
		SPLOG("Flyweight validation complete");
	}
}

manifest mem_seqset_task::build_seqset()
{
	manifest seqset_manifest;

	subprogress sp_load([this](double p) { update_progress(p); }, 0.0, 0.1);
	subprogress sp_compute([this](double p) { update_progress(p); }, 0.1, 1.0);
	SPLOG("Running with %lu threads, max_buf_size: %lu", num_threads, m_max_buf_size);

	reference ref(ref_name);
	load_repo(ref, sp_load);
	seqset_manifest = do_mem_seqset(sp_compute);

	return seqset_manifest;
}

size_t mem_seqset_task::count_flyweights() const
{
	size_t flywt_count = 0ULL;

	if (is_paired) {
		flywt_count += input.get_num_records() * 4; // Forward, reverse, mate forward, mate reverse.
	} else {
		flywt_count += input.get_num_records() * 2; // Forward, reverse.  No mate pairs.
	}

	SPLOG("mem_seqset_task::count_flyweights> %lu", flywt_count);
	return flywt_count;
}

void mem_seqset_task::run()
{
	validate();
	set_output(build_seqset());
}

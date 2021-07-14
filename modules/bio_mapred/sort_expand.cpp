#include "modules/bio_base/seqset.h"
#include "modules/mapred/sort_task.h"
#include "modules/bio_mapred/sort_expand.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/io/parallel.h"

#include <boost/bind.hpp>
#include <boost/dynamic_bitset.hpp>
#include <parallel/algorithm>

REGISTER_TASK(sort_expand_task);
REGISTER_TASK(sort_expand_part_task);

static dna_slice decode_flyweight(const dna_sequence& repo, uint32_t read_size, uint32_t f) {
	uint32_t read = ((f/2) / read_size);
	uint32_t size = read_size - ((f/2) % read_size); 
	auto start = ((f%2==0) ? 
			repo.begin() + (read * read_size + (read_size - size)) :
			(repo.begin() + (read * read_size + size - 1)).rev_comp());
	return dna_slice(start, size_t(size));
}

void sort_expand_part_task::run() 
{
	progress_handler_t my_prog = boost::bind(&sort_expand_part_task::void_progress, this, _1);
	subprogress sp1(my_prog, 0.0, 0.4); // Pull in actual reads
	subprogress sp2(my_prog, 0.7, 1.0); // Write out

	std::string key;
	corrected_reads value;
	size_t num_records = input.num_records;
	
	SPLOG("Reading initial read");	
	auto source = input.build();
	if (!source->read_msgpack(key, value)) {
		throw io_exception("Trying to make empty PBWT");
	}
	size_t read_size = value[0].corrected.size();
	size_t num_seqs = num_records * read_size * 2;
	size_t c = 0;
	SPLOG("Read size = %zu, number of reads = %zu, reading them in", read_size, num_records);	

	dna_sequence repo(read_size * num_records);
	std::vector<uint32_t> flyweights;
	flyweights.reserve(num_seqs);
	do {
		for(size_t i = 0; i < read_size; i++) {
			repo[c*read_size + i] = value[0].corrected[i];
			flyweights.push_back(2*(c*read_size + i));
			flyweights.push_back(2*(c*read_size + i) + 1);
		}
		c++;
		sp1(double(c) / double(num_records));
	} while(source->read_msgpack(key, value));

	SPLOG("Doing the (parallel?) sort");
	
        auto sort_future = std::async(std::launch::async, [&] {
            parallel_sort_in_place(
			flyweights.begin(), flyweights.end(), [&](const uint32_t &a, const uint32_t &b) {
			dna_slice as = decode_flyweight(repo, read_size, a);
			dna_slice bs = decode_flyweight(repo, read_size, b);
			return as < bs;
		});
        });
	future_watchdog(sort_future, my_prog, .5);

	SPLOG("Doing the (parallel?) compare");
	uint32_t *buf = new uint32_t[(num_seqs + 31) / 32];
	memset(buf, 0, sizeof(uint32_t) * ((num_seqs + 31) / 32));
	
        auto compare_future = std::async(std::launch::async, [&] {
		#ifdef _OPENMP	
		#pragma omp parallel for
		#endif
		for(size_t i = 0; i < num_seqs - 1; i++) {
			dna_slice as = decode_flyweight(repo, read_size, flyweights[i]);
			dna_slice bs = decode_flyweight(repo, read_size, flyweights[i+1]);
			if (bs.size() < as.size() || bs.subseq(0, as.size()) != as) {
				__sync_fetch_and_or(&buf[i / 32], 1 << (i % 32)); 
			}
		}
		buf[(num_seqs -1) / 32] |= (1 << ((num_seqs - 1) % 32));
	});
	future_watchdog(compare_future, my_prog, .6);

	SPLOG("Writing data");
	manifest mout;
	auto out = output.build(get_root(), "exps",  mout);
	size_t tot = 0;
	for(size_t i = 0; i < num_seqs; i++) {
		sp2(double(i) / double(num_seqs));
		if (buf[i / 32] & (1 << (i % 32)))
		{
			tot++;
			dna_slice as = decode_flyweight(repo, read_size, flyweights[i]);
			dna_sequence seq(as.begin(), as.end());  // TODO: Remove needless copy
			//SPLOG("Kept: %s", seq.as_string().c_str());
			out->write_msgpack(seq, 0);
		} else {
			//dna_slice as = decode_flyweight(repo, read_size, flyweights[i]);
			//SPLOG("Dropped: %s", as.as_string().c_str());
		}
	}
	SPLOG("Kept %zu sequences", tot);

	delete[] buf;
	out->close();
	set_output(mout);
}

void sort_expand_task::run() {
	if (m_state == 0) {
		split_progress(0.01, 0.25);
		std::vector<input_stream_params> inputs;
		input.split_by_goal_size(inputs, split_size);

		for (size_t i = 0; i < inputs.size(); i++) {
			std::unique_ptr<sort_expand_part_task> t = make_unique<sort_expand_part_task>();
			t->input = inputs[i];
			t->output.presorted = true;
			t->output.sort = "prefix";
			t->output.encoding = "null";
			t->output.goal_size = 256*1024*1024;
			m_subtasks.push_back(add_subtask(std::move(t)));
			update_progress(double(i) / double(inputs.size()));
		}
		m_state = 1;
	} else if (m_state == 1) {
		split_progress(0.01, 0.01);
		manifest out("prefix");

		for (size_t i = 0; i < m_subtasks.size(); i++) {
			manifest subout;
			get_output(subout, m_subtasks[i]);
			out.add(subout);
			update_progress(double(i) / double(m_subtasks.size()));
		}
		m_subtasks.clear();
		std::unique_ptr<sort_task> t = make_unique<sort_task>();
		t->input = out;
		t->is_summary = true;
		t->reduce = "prefix";
		m_subtasks.push_back(add_subtask(std::move(t)));
		m_state = 2;
	} else if (m_state == 2) {
		manifest out;
		get_output(out, m_subtasks[0]);
		// Update sort to prevent further prefix merging
		out.set_sort("dna");
		set_output(out);
	}
}

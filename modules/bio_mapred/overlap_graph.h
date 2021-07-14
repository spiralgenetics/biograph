
#pragma once

#include <boost/format.hpp>
#include "modules/bio_base/dna_sequence.h"
#include <algorithm>
#include <vector>
#include "modules/io/log.h"
#include "modules/io/track_mem.h"

/* Concepts used by read overlap graph

// A vector like thing that knows its size and had dna_slices in it
class read_vector {
	size_t size() const;
	const dna_slice& operator[](size_t i) const;
}

// A progess handler, called with x = 0.0 at start to x = 1.0 at end.
void progress(double x);

// A function like callback object that is called to alert the caller that
// two 'reads' are duplicates, and read 'i' will be kept, and read 'j' (the dup)
// will be removed from the list.  Flipped indicate if j is a reverse complement dup
// or a straight dup.
void on_duplicate(size_t i, size_t j, bool flipped);

// A function like callback object that reports 'overlap' results.
// i is the read index that overlaps, flipped is true if the overlap is for the reverse
// complement of read i, and overlap is the number of matching bases.
void on_overlap(size_t i, bool flipped, uint8_t overlap)

*/


template<class ReadVec>
class overlap_graph
{
public:
	typedef ReadVec read_vec_t;

private:
	struct read_index {
		read_index() {}
		read_index(size_t _index, bool _flipped)
			: index(_index), flipped(_flipped) {}
		unsigned int index : 31;
		unsigned int flipped : 1;
	};
	struct seq_range {
		seq_range(const dna_const_iterator& _it, size_t _len)
			: it(_it), len(_len)
		{}
		seq_range right(size_t offset) { return seq_range(it + offset, len - offset); }
		dna_const_iterator it;
		size_t len;
	};
	const read_vec_t& m_reads;
	tracked_vector<read_index> m_lookup;
	bool priv_less_than_inc_index(const read_index& i, const read_index& j) {
		dna_slice si = m_reads[i.index];
		dna_slice sj = m_reads[j.index];
		dna_const_iterator iti = (i.flipped ? si.rcbegin() : si.begin());
		dna_const_iterator itj = (j.flipped ? sj.rcbegin() : sj.begin());
		dna_compare_result r = subseq_compare(iti, itj, si.size(), sj.size());
		if (r != dna_compare_result::EQUAL) return r == dna_compare_result::FIRST_IS_LESS || r == dna_compare_result::FIRST_IS_PREFIX;
		if (i.index != j.index) return i.index < j.index;
		return i.flipped < j.flipped;
	}
	bool priv_less_than(const read_index& i, const read_index& j) {
		dna_slice si = m_reads[i.index];
		dna_slice sj = m_reads[j.index];
		dna_const_iterator iti = (i.flipped ? si.rcbegin() : si.begin());
		dna_const_iterator itj = (j.flipped ? sj.rcbegin() : sj.begin());
		return subseq_lessthan(iti, itj, si.size(), sj.size());
	}
	bool priv_less_than(const read_index& i, const seq_range& j) {
		dna_slice si = m_reads[i.index];
		dna_const_iterator iti = (i.flipped ? si.rcbegin() : si.begin());
		return subseq_lessthan(iti, j.it, si.size(), j.len);
	}
	bool priv_front_match(const read_index& i, const seq_range& j) {
		dna_slice si = m_reads[i.index];
		dna_const_iterator iti = (i.flipped ? si.rcbegin() : si.begin());
		return subseq_equal(iti, j.it, j.len);
	}
	std::string nice_seq(const read_index& i) {
		dna_slice s = m_reads[i.index];
		if (i.flipped) s = s.rev_comp();
		return printstring("%d:%d:%s", i.index, i.flipped, s.as_string().c_str());
	}
public:
	// Constructs an overlap graph
  overlap_graph(const read_vec_t& reads) : m_reads(reads), m_lookup(track_alloc("overlap_graph:lookup")) {}

	// Prepares an overlap graph for use, reports progress and finds duplicates
	template<typename ProgressType, typename OnDuplicateType>
	void prepare(ProgressType progress, OnDuplicateType on_duplicate)
	{
		// Setup some convenience variables
		size_t size = m_reads.size();
		double dsize = double(2*size);
		auto less_lambda = [this](const read_index& i, const read_index& j) {
			return this->priv_less_than_inc_index(i, j);
		};
		// Early exit for pathological case
		if (size == 0) return;

		// Allocate some space
		m_lookup.resize(2 * size);

		// First we "make" a heap, one element always is a heap.
		m_lookup[0] = read_index(0, 0);

		// Do a 'heap_sort' while tracking progress
		// Push items onto a heap
		SPLOG_P(LOG_DEBUG, "overlap_graph::prepare> Pushing into the heap");
		for(size_t i = 1; i < 2*size; i++) {
			progress(.1 * double(i) / dsize);
			m_lookup[i] = read_index(i / 2, (i % 2 == 1));
			push_heap(m_lookup.begin(), m_lookup.begin() + i, less_lambda);
		}

		// Pop items from heap
		SPLOG_P(LOG_DEBUG, "overlap_graph::prepare> Popping from the heap");
		for(size_t i = 0; i < 2*size; i++) {
			progress(.1 + .8 * double(i) / dsize);
			pop_heap(m_lookup.begin(), m_lookup.begin() + (2*size - i), less_lambda);
		}

		// Do in place dedup
		SPLOG_P(LOG_DEBUG, "overlap_graph::prepare> Doing dedup");
		size_t start = 0;  // Start of current 'group'
		size_t out = 1;	 // Place to write next output,	start is already written
		size_t dup_count = 0;
		for(size_t i = 1; i < 2*size; i++) {
			progress(.9 + .1 * double(i) / dsize);
			read_index& sri = m_lookup[start];
			read_index& iri = m_lookup[i];
			if (priv_less_than(sri, iri)) {
				// Non dup, 'copy' across
				m_lookup[out++] = iri;
				start = i;
			} else {
				// Dup, alert only once, and make sure to skip
				// call on 'palendromic' reads
				if (sri.flipped == 0 && sri.index != iri.index)
					on_duplicate(sri.index, iri.index, iri.flipped);
				dup_count++;
			}
		}
		SPLOG("%s", str(boost::format("overlap_graph::prepare> %1% reads passed and %2% reads marked as duplicates.")
			% out % dup_count).c_str());

		// Remove extra space from m_lookup
		m_lookup.resize(out);
		m_lookup.shrink_to_fit();
	}

	// Given a read index, i, finds all overlaps reads in the forward (fwd) or reverse
	// direction, which overlap by at least min_overlap, and calls the on_overlap
	// function object to report them.  Always returns overlaps in the order of most
	// overlap to least overlap
	template<typename OverlapType>
	void find_overlaps(size_t i, bool fwd, size_t min_overlap, OverlapType on_overlap)
	{
		auto less_lambda = [this](const read_index& i, const seq_range& j) {
			return this->priv_less_than(i, j);
		};
		seq_range sr(fwd ? m_reads[i].begin() : m_reads[i].rcbegin(), m_reads[i].size());
		for(uint8_t overlap = sr.len - 1; overlap >= min_overlap; overlap--) {
			seq_range sr2 = sr.right(sr.len - overlap);
			auto it = lower_bound(m_lookup.begin(), m_lookup.end(), sr2, less_lambda);
			while((it != m_lookup.end()) && priv_front_match(*it, sr2)) {
				on_overlap(it->index, (fwd ? it->flipped : !it->flipped), overlap);
				it++;
			}
		}
	}
};

#pragma once

#include <map>
#include <vector>
#include "modules/bio_base/dna_sequence.h"

class ipileup {
 public:
  virtual ~ipileup() = default;
  virtual size_t depth_at(size_t position) const = 0;
  virtual size_t fwd_at(size_t position) const = 0;
  virtual size_t tot_qual_at(size_t position) const = 0;
};

// The pileup class represents a way to collect information about
// a set of reads and how they relate to a given sequence.  We align each
// read to its position on the sequence and then collect statistics
// When aligning the read, we choose postion to minimize phred error score
// of mismatched bases, and don't current support inserts or deletes,
// which aren't needed for anchored assembly
class pileup : public ipileup
{
	struct base_info {
		uint32_t count = 0;
		uint32_t tot_qual = 0;
		uint32_t fwd = 0;
	};
	typedef std::map<dna_base, base_info> base_quals_t;
	typedef std::vector<base_quals_t> pileup_t;
	struct read_info {
		size_t offset;
		std::string name;
		dna_sequence seq;
		std::string qual;
	};
public:
	// Construct a new and intially empty pileup.
	pileup(const dna_sequence& sequence, // The sequence to align reads to
		size_t max_cost);            // The maximum cost for read alignment

	// Returns the sequence a pileup is piled up against
	const dna_sequence& get_sequence() { return m_sequence; }

	// Add a read to the pileup, returns pos or -1 if cost is too high
	int add_read(const std::string& name, // The name of the read
		const dna_sequence& read_seq,  // The read sequence
		const std::string quality, // The quality (as phred scores)
		bool fwd,  // Is the reads as is, or reversed
		int pos);  // The position (if known)

	// Return depth of matching reads at position
	size_t depth_at(size_t position) const override;
	size_t fwd_at(size_t position) const override;
	size_t tot_qual_at(size_t position) const override;

	// Finds all reads that overlap a given region
	// Start and end are given in C++ iterator style [start, end)
	// as offsets from the start of sequence
	// std::vector<read_info> relevant_reads(size_t start, size_t end);

	void print();
private:
	int best_match_pos(const dna_sequence& read_seq, const std::string quality);

	dna_sequence m_sequence;
	size_t m_max_cost;
	pileup_t m_pileup;		
	//std::vector<read_info> m_reads;		
};



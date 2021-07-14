#pragma once

#include "modules/io/transfer_object.h"
#include "modules/io/utils.h"
#include "modules/bio_base/seq_position.h"

// Let's talk terms and coordinate systems.
//
// FASTA files contain one reference assembly.
//
// The reference assembly is a collection of one or more scaffolds.  Each
// scaffold contains one or more supercontigs (contiguous bases of ACTG)
// likely separated by regions of N bases:
//
// >1 | The first scaffold (sometimes it's a chromosome) | more junk here
// NNNACTGACTGNNNNCGGTA
// CCGNNNATTACANGGATNNN
// >2 | The next scaffold
// NNNNNNCCATGCCANNNNTT
// GTTACCCATGNNNNCTTTAT
//
// The scaffold name is everything after the > up to the first whitespace
// character. Everything from the first whitespace forward is a comment and is
// ignored.
//
// There are 40 nucleotide bases in the first scaffold. If the first base in
// this scaffold is at index 0, each base would have a position of:
//
// N  N  N  A  C  T  G  A  C  T  G  N  N  N  N  C  G  G  T  A
// 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
// C  C  G  N  N  N  A  T  T  A  C  A  N  G  G  A  T  N  N  N
// 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
//
// There are four supercontigs in the sequence, at 1:3, 1:15,
// 1:26, and 1:33.
//
// The indices in this coordinate system are referred to as the sequence
// position. They are expressed as seq_position objects, with a scaffold id
// and a position (see seq_position.h)
//
// If we remove the N bases, the indices look like this:
//
//          A  C  T  G  A  C  T  G              C  G  G  T  A
//          3  4  5  6  7  8  9  10             15 16 17 18 19
// C  C  G           A  T  T  A  C  A     G  G  A  T
// 20 21 22          26 27 28 29 30 31    33 34 35 36
//
// Grouping the sequences by supercontig, we have:
//
// 1:3
// A  C  T  G  A  C  T  G
// 3  4  5  6  7  8  9  10
// 1:15
// C  G  G  T  A  C  C  G
// 15 16 17 18 19 20 21 22
// 1:26
// A  T  T  A  C  A
// 26 27 28 29 30 31
// 1:33
// G  G  A  T
//
// If we indexed all of those positions starting at 0, it would look like this:
//
//    S/C:  1:3                     1:15                    1:26              1:33
//   BASE:  A  C  T  G  A  C  T  G  C  G  G  T  A  C  C  G  A  T  T  A  C  A  G  G  A  T
// SEQPOS:  3  4  5  6  7  8  9  10 15 16 17 18 19 20 21 22 26 27 28 29 30 31 33 34 35 36
//  INDEX:  0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
//
// The indices in this discontiguous coordinate system are called the flat
// position. For convenience, we process DNA sequences internally using 0-indexed
// flat space.
//
// To translate from flat position to sequence position, use reference->get_seq_position()
// (which takes a size_t flat position and returns a seq_position).
//
// To translate from a sequence position to flat position, use reference->flatten()
// (which takes a seq_position or (scaffold name, offset) and returns a size_t).
//
// Note that humans expect bases to be 1-indexed, in seq_position space. But
// we work internally in 0-indexed, flat position space. This is fixed up in
// the API by using get_seq_position() and subtracting or adding 1 when
// needed.
//

struct supercontig
{
	supercontig() : offset(0), len(0) {} // For deserialization
	supercontig(const std::string& _scaffold, size_t _offset, size_t _len);
	TRANSFER_OBJECT { VERSION(0); FIELD(scaffold_name); FIELD(name); FIELD(offset); FIELD(len); }
	bool operator<(const supercontig& rhs) const;

	std::string scaffold_name;
	std::string name;
	size_t offset;   // Offset in scaffold
	mutable size_t tot_offset;  // Offset in global flattened order
	size_t len;
};

struct scaffold
{
	scaffold() : len(0) {} // For deserialization
	scaffold(const std::string& _name, size_t _len, int index);
	TRANSFER_OBJECT { VERSION(0); FIELD(name); FIELD(len); FIELD(index);}
	bool operator<(const scaffold& rhs) const;
	std::string name;
	size_t len;
	int index;
	long start;
};

struct reference_assembly
{
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(supercontigs);
		FIELD(scaffolds);
		FIELD(scaffold_order);
		if (IS_DESERIALIZE)
			generate_tables();
	}

	std::set<supercontig> supercontigs;
	std::set<scaffold> scaffolds;
	std::vector<std::string> scaffold_order;

	size_t size() const;

	const scaffold& get_scaffold(const std::string& name) const;
	const supercontig& get_supercontig(const std::string& name) const;
	const std::vector<std::string>& get_supercontig_order() const;
	seq_position get_seq_position(size_t pos) const;
	const supercontig& get_supercontig(size_t pos) const;

	size_t flatten(const seq_position& loc) const;
	size_t flatten(std::string scaffold_name, size_t pos) const;
	std::pair<size_t, size_t> flatten_range(
		const std::string& contig_name
		, unsigned long start
		, unsigned long end
		, bool use_exact_loci = true
	) const;

private:
	mutable size_t m_size;
	mutable std::vector<std::string> m_supercontig_order;
	mutable std::vector<size_t> m_supercontig_start;
	mutable std::map<std::string, uint32_t> m_supercontig_offset;

public:
	// The supercontig order and base convert tables
	// are not stored in JSON, but is built during deserialize
	void generate_tables() const;
};

struct flatten_exception : public io_exception
{
	flatten_exception(const std::string& message) : io_exception(message) {}
	flatten_exception(const boost::format& formatted_message) : io_exception(formatted_message.str()) {}
};

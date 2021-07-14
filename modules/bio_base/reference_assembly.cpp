#include <limits>
#include <algorithm>

#include "modules/io/log.h"
#include "modules/bio_base/reference_assembly.h"

// A supercontig is a contiguous region of nucleotides (non-N regions of a FASTA)
supercontig::supercontig(const std::string& _scaffold, size_t _offset, size_t _len)
	: scaffold_name(_scaffold)
	, name(printstring("%s:%d", _scaffold.c_str(), (int) _offset))
	, offset(_offset)
	, len(_len)
{}

bool supercontig::operator<(const supercontig& rhs) const
{
	if (scaffold_name != rhs.scaffold_name)
		return scaffold_name < rhs.scaffold_name;
	return offset < rhs.offset;
}

// A scaffold is a region of supercontigs, separated by N regions.
scaffold::scaffold(const std::string& _name, size_t _len, int _index)
	: name(_name)
	, len(_len)
	, index(_index)
{}

bool scaffold::operator<(const scaffold& rhs) const
{
	return name < rhs.name;
}

// A reference_assembly is a group of scaffolds representing a reference organism.
size_t reference_assembly::size() const
{
	return m_size;
}

// Return a scaffold object with the given name. Scaffold names cannot contain whitespace.
const scaffold& reference_assembly::get_scaffold(const std::string& name) const
{
	scaffold key(name, 0, 0);
	std::set<scaffold>::const_iterator it = scaffolds.find(key);
	if (it == scaffolds.end())
		throw io_exception(printstring("Lookup of invalid scaffold %s", name.c_str()));
	return *it;
}

// NOTE: reference_assembly::get_supercontig() is overloaded.
// Return the supercontig containing the given position (in flattened coordinate space)
const supercontig& reference_assembly::get_supercontig(size_t pos) const
{
	auto it = upper_bound(m_supercontig_start.begin(), m_supercontig_start.end(), pos);
	it--;
	size_t index = it - m_supercontig_start.begin();
	const std::string& name = m_supercontig_order[index];
	return get_supercontig(name);
}

// Return the supercontig object with the given name "scaffold:start_pos" (eg. "1:227417")
const supercontig& reference_assembly::get_supercontig(const std::string& name) const
{
	std::string::size_type it = name.rfind(':');
	if (it == std::string::npos)
		throw io_exception("get_supercontig: Missing : in name '" + name + "'");
	std::string scaf = name.substr(0,it);

	size_t sc_offset = atoi(name.substr(it + 1).c_str());
	supercontig key(scaf, sc_offset, 0);

	std::set<supercontig>::const_iterator it2 = supercontigs.find(key);
	if (it2 == supercontigs.end())
		throw io_exception("get_supercontig: Cannot find supercontig '" + name + "'");
	return *it2;
}

// Natural order of all supercontigs.
const std::vector<std::string>& reference_assembly::get_supercontig_order() const
{
	return m_supercontig_order;
}

// Given a flattened position, return a seq_position (unflatten).
seq_position reference_assembly::get_seq_position(size_t pos) const
{
	const supercontig& sc = get_supercontig(pos);
	long seq_pos = pos - sc.tot_offset;
	int scaffold_id = get_scaffold(sc.scaffold_name).index;
	seq_pos += sc.offset;
	seq_position sp(scaffold_id, seq_pos);
	return sp;
}

// Given a sequence position, return the flattened position
//
// Flattened coordinate space is 0-based, discontiguous, and only includes DNA
// regions (no N blocks). See reference_assembly.h.
//
size_t reference_assembly::flatten(const seq_position& loc) const
{
	std::string scaffold_name = scaffold_order[loc.scaffold_id];
	long pos = loc.position;

	if (supercontigs.size() == 0) { // Lots of edge cases, early exit to avoid
		throw io_exception("The reference does not appear to contain any sequence data");
	}

	supercontig key(scaffold_name, pos, 0);
	// Find first element > key
	std::set<supercontig>::const_iterator it = supercontigs.upper_bound(key);
	if (it == supercontigs.begin()) { // Before first element, forget it
		throw flatten_exception(boost::format(
			"%1%:%2% comes before the first non-N sequence in the reference (offset %3%)"
			) % scaffold_name % pos % it->offset);
	}

	--it;  // Find last element <= key
	if ((long) it->offset > pos) {
		throw flatten_exception(boost::format(
			"%1%:%2% is in a part of the reference (offset %3%) with no sequence data (only N bases)"
			) % scaffold_name % pos % it->offset);
	}

	if (pos - it->offset > it->len) {
		throw flatten_exception(boost::format(
				"%1%:%2% is in a part of the reference (offset %3%, length %4%) with no sequence data (only N bases)"
				) % scaffold_name % pos % it->offset % it->len);
	}

	// Correct coordinates
	pos -= it->offset;

	return it->tot_offset + size_t(pos);
}

// Given a scaffold name and offset, return the flattened position
size_t reference_assembly::flatten(std::string scaffold_name, size_t pos) const
{
	return flatten(seq_position(get_scaffold(scaffold_name).index, pos));

}

// Populate the m_supercontig_* lists.
// These should really be called m_contig_*
void reference_assembly::generate_tables() const
{
	uint32_t offset = 0;
	for (std::string scaffold_name : scaffold_order)
	{
		auto starting_supercontig_it = supercontigs.lower_bound(supercontig(scaffold_name, 0, 0));
		auto ending_supercontig_it = supercontigs.upper_bound(supercontig(scaffold_name, std::numeric_limits<std::size_t>::max(), 0));

		for(auto it = starting_supercontig_it; it != ending_supercontig_it; ++it)
		{
			m_supercontig_order.push_back(it->name);
			m_supercontig_start.push_back(offset);
			it->tot_offset = offset;
			offset += it->len;
		}
	}
	m_size = offset;
}

// Given a contig name, start, and end (in absolute coordinates), return a pair of <start, end>
// in flat space.
//
// Throws if start and end are not a valid address, or if they cross supercontig boundaries.
//
// If use_exact_loci is false, attempt to move start or end (not both) outside of an N region.
//
std::pair<size_t, size_t> reference_assembly::flatten_range(
	const std::string& contig_name
	, unsigned long start
	, unsigned long end
	, bool use_exact_loci
) const
{
	int contig_id = get_scaffold(contig_name).index;

	std::string error_message;
	const size_t k_bad_locus = std::numeric_limits<size_t>::max();
	size_t flat_start = k_bad_locus;
	size_t flat_end = k_bad_locus;

	try {
		flat_start = flatten(seq_position(contig_id, start));
	} catch (const flatten_exception& e) {
		if (use_exact_loci) {
			throw;
		}
		error_message = e.what();
	}

	try {
		flat_end = flatten(seq_position(contig_id, end));
	} catch (const flatten_exception& e) {
		if (use_exact_loci) {
			throw;
		}

		if (flat_start == k_bad_locus) {
			throw io_exception(boost::format(
				"Both the start (%1%:%2%) and end (%1%:%3%) loci are in reference areas"
				" with no data (N bases).\n\nStart: %4%\n\nEnd: %5%")
				% contig_name % start % end % error_message % e.what());
		}
	}

	if (start >= end) {
		throw io_exception(boost::format(
			"The start location \"%1%:%2%\" must come before the end location \"%1%:%3%\""
			" on the reference.") % contig_name % start % end);
	}

	if (flat_start == k_bad_locus) {	// Start is bad, but end is good
		const supercontig& end_supercontig = get_supercontig(flat_end - 1);
		flat_start = flatten(seq_position(contig_id, end_supercontig.offset));
		SPLOG("reference_wrapper::make_range> start locus is not in supercontig."
			"  Adjusting start forward to beginning of supercontig at %s:%lu"
			, contig_name.c_str(), get_seq_position(flat_start).position);
	} else if (flat_end == k_bad_locus) {	// End is bad, but start is good
		const supercontig& start_supercontig = get_supercontig(flat_start);
		flat_end = flatten(seq_position(contig_id, start_supercontig.offset + start_supercontig.len));
		SPLOG("reference_wrapper::make_range> end locus is not in supercontig."
			"  Adjusting end back to end of supercontig at %s:%lu"
			, contig_name.c_str(), get_seq_position(flat_end - 1).position + 1);
	}

	if (get_supercontig(flat_start).name != get_supercontig(flat_end - 1).name) {
		throw io_exception(boost::format(
			"The start location \"%1%:%2%\" and end location \"%1%:%3%\""
			" are not in the same contiguous region. (%4% vs. %5%)")
			% contig_name % start % end % get_supercontig(flat_start).name % get_supercontig(flat_end - 1).name);
	}

	return std::make_pair(flat_start, flat_end);
}

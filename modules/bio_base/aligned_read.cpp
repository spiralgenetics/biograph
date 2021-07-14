#include "base/base.h"
#include "modules/bio_base/aligned_read.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

namespace
{
	void parse_sam_optional_field(aligned_read& the_read, const std::string& the_field);
}

// SAM fields are documented here: https://samtools.github.io/hts-specs/SAMv1.pdf
//
// Field names are listed below with zero indexing.
//
//  0: QNAME
//  1: FLAG
//  2: RNAME
//  3: POS
//  4: MAPQ
//  5: CIGAR
//  6: RNEXT
//  7: PNEXT
//  8: TLEN
//  9: SEQ
// 10: QUAL

bool parse_sam(const reference_assembly& ref, aligned_read& out, const std::string& sam_line)
{
	std::vector<std::string> fields;
	boost::split(fields, sam_line, boost::is_any_of("\t"));
	if (fields.size() < 11) {
		return false;
	}
	out.read_name = fields[0];
	out.flags = atoi(fields[1].c_str());
	if (fields[6] == "=") {
		fields[6] = fields[2];
	}
	if (fields[2] == "*" || atoi(fields[3].c_str()) == 0) {
		out.ref_pos = seq_position();
	}
	else {
		long pos = atoi(fields[3].c_str());
		out.ref_pos = ref.get_seq_position(ref.flatten(fields[2], pos));
	}
	out.map_quality = atoi(fields[4].c_str());
	out.cigar = fields[5];
	if (fields[6] == "*" || atoi(fields[7].c_str()) == 0) {
		out.mate_pos = seq_position();
	}
	else {
		long pos = atoi(fields[7].c_str());
		out.mate_pos = ref.get_seq_position(ref.flatten(fields[6], pos));
	}
	out.tlen = atoi(fields[8].c_str());
	out.seq = fields[9];
	out.qual = fields[10];

	for(unsigned int i = 11; i < fields.size(); i++) {
		parse_sam_optional_field(out, fields[i]);
	}

	return true;
}

std::string print_sam(const reference_assembly& ref, const aligned_read& in, bool use_supercontig_coords)
{
	std::string ref_name = "*";
	std::string mate_name = "*";
	long ref_pos = 0;
	long mate_pos = 0;

	if(use_supercontig_coords) {
		throw io_exception("Supercontig coordinates are unsupported.");
	}
	if (in.ref_pos.valid()) {
		ref_name = ref.get_supercontig_order()[in.ref_pos.scaffold_id];
		ref_pos = in.ref_pos.position;
	}
	if (in.mate_pos.valid()) {
		mate_name = ref.get_supercontig_order()[in.mate_pos.scaffold_id];
		mate_pos = in.mate_pos.position;
		if (mate_name == ref_name)
			mate_name = "=";
	}
	return printstring(
		"%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%ld\t%s\t%s\t%s",
		in.read_name.c_str(),
		in.flags,
		ref_name.c_str(),
		ref_pos,
		in.map_quality,
		in.cigar.c_str(),
		mate_name.c_str(),
		mate_pos,
		in.tlen,
		in.seq.c_str(),
		in.qual.c_str(),
		"RG:Z:Spiral" // HACK:  Print a fake read group to make GATK happy.
	);

}


namespace
{
	void parse_sam_optional_field(aligned_read& the_read, const std::string& the_field)
	{
		static const boost::regex field_regex("([A-Za-z][A-Za-z0-9]):([AifZHB]):(.+)");
		boost::smatch the_match;

		if (! boost::regex_match(the_field, the_match, field_regex)) {
			std::string error_string("A sam optional field failed to parse: '");
			error_string += the_field;
			error_string += "' on read ID '";
			error_string += the_read.read_name;
			error_string += "'";
			throw io_exception(error_string.c_str());
		}

		std::string tag(the_match[1].first, the_match[1].second);
		std::string type(the_match[2].first, the_match[2].second);
		std::string value(the_match[3].first, the_match[3].second);

		if (tag == "RG") {
			CHECK_EQ(type, "Z");
			the_read.read_group_id = value;
		}
	}
} //	End blank namespace

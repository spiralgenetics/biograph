
#ifndef __unaligned_read_h__
#define __unaligned_read_h__

#include <boost/container/small_vector.hpp>

#include "modules/io/transfer_object.h"
#include "modules/io/keyvalue.h"
#include "modules/bio_base/seq_position.h"

struct read_id
{
	static const int type_id = 'u';
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(pair_name);
	}
	std::string pair_name;
};

inline bool operator!=(const read_id& lhs, const read_id& rhs)
{
	return lhs.pair_name != rhs.pair_name;
}

struct unaligned_read
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(pair_number);
		FIELD(name_suffix);
		FIELD(sequence);
		FIELD(quality);
		FIELD(original_sequence_id);
		FIELD(ref_loc);
		FIELD(mismatches, -1);
	};
	int pair_number;
	std::string name_suffix;
	std::string original_sequence_id;
	std::string sequence;
	std::string quality;	
	seq_position ref_loc;
	int mismatches = -1;
	void convert_phred64();
	void trim3(int len);
	void trim5(int len);
};

using unaligned_reads = boost::container::small_vector<unaligned_read, 2 /* 2 mates per pair */>;
SET_TYPE_ID(unaligned_reads, 'U')

void parse_read_name(const std::string& name, std::string& pair_name, unaligned_read& out);
std::string build_read_name(const std::string& pair_name, const unaligned_read& out);

void convert_phred64(unaligned_reads& reads);
void trim3(unaligned_reads& reads, int len);
void trim5(unaligned_reads& reads, int len);

#endif


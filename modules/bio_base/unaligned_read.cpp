
#include "modules/bio_base/unaligned_read.h"
#include "modules/io/utils.h"
#include <boost/regex.hpp>

void unaligned_read::convert_phred64()
{
	for(size_t i = 0; i < quality.size(); i++)
	{
		quality[i] -= 31;
		if (quality[i] <= 32)
			throw io_exception("phred_64 option set and quality values < 64");
	}
}

void unaligned_read::trim3(int trim)
{
	int new_size = std::max(0, (int) sequence.size() - trim);
	sequence = sequence.substr(0, new_size);
	quality = quality.substr(0, new_size);
}

void unaligned_read::trim5(int trim)
{
	int remove = std::min(trim, (int) sequence.size());
	sequence = sequence.substr(remove);
	quality = quality.substr(remove);
}

void convert_phred64(unaligned_reads& reads)
{
	for(size_t i = 0; i < reads.size(); i++)
		reads[i].convert_phred64();
}

void trim3(unaligned_reads& reads, int trim)
{
	for(size_t i = 0; i < reads.size(); i++)
		reads[i].trim3(trim);
}

void trim5(unaligned_reads& reads, int trim)
{
	for(size_t i = 0; i < reads.size(); i++)
		reads[i].trim5(trim);
}

void parse_read_name(const std::string& name, std::string& key, unaligned_read& read)
{
	key = name;
	read.pair_number = 0;
	read.name_suffix = "";
	read.original_sequence_id = name;
}

std::string build_read_name(const std::string& key, const unaligned_read& read)
{
	if(read.original_sequence_id.empty())
	{
		if (read.pair_number == 0)
			return key;
		else
			return printstring("%s%d%s", key.c_str(), read.pair_number, read.name_suffix.c_str());
	}
	else
	{// new code for the new uploaded paired data.
		return read.original_sequence_id;
	}
}

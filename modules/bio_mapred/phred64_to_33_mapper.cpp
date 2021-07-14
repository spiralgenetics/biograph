#include "modules/io/msgpack_transfer.h"
#include "modules/io/log.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_mapred/phred64_to_33_mapper.h"

REGISTER_1(mapper, phred64_to_33, const std::string&);

class phred_exception
{
public:
	phred_exception(char bad_quality) : m_bad_quality_char(bad_quality) {}
	char get_bad_quality() const { return m_bad_quality_char; }
private:
	char m_bad_quality_char;
};
        
void phred64_to_33_mapper::map(const std::string& key, const std::string& value, kv_sink& context)
{
	//SPLOG("phred64_to_33_mapper::map got key %s", key.c_str());
	read_id the_read_id;
	unaligned_reads the_reads;
	
	try
	{
		msgpack_deserialize(the_read_id, key);
		msgpack_deserialize(the_reads, value);
	}
	catch(const deserialization_error& e)
	{
		std::string error_string("in phred64_to_33_mapper: ");
		error_string += e.message();
		error_string += ", ";
		error_string + the_read_id.pair_name;
		throw io_exception(error_string); 
	}
	
	std::for_each(the_reads.begin(), the_reads.end(),
		[this](unaligned_read& a_read) mutable
		{
			try
			{
				convert_64_to_33(a_read.quality);
			}
			catch (phred_exception& pe)
			{
				throw io_exception(make_exception_string(a_read, pe.get_bad_quality()));
			}
		}
	);
	
	context.write(msgpack_serialize(the_read_id), msgpack_serialize(the_reads));
}

void phred64_to_33_mapper::convert_64_to_33(std::string& qualities_string) const
{
	std::for_each(qualities_string.begin(), qualities_string.end(),
		[](char& a_quality) mutable
		{
			if (a_quality < 64) throw phred_exception(a_quality);
			a_quality -= (64 - 33);
		}
	);
}

std::string phred64_to_33_mapper::make_exception_string(const unaligned_read& a_read, char bad_quality) const
{
	std::string error_string("in phred64_to_33_mapper: Read ID ");
	error_string += a_read.original_sequence_id;
	error_string += " has quality \"";
	error_string += bad_quality;
	error_string += "\" which is illegal in phred-64.";
	
	return error_string;
}

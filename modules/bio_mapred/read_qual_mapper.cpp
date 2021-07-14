
#include "modules/bio_mapred/read_qual_mapper.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/io/msgpack_transfer.h"

REGISTER_1(mapper, read_qual, const std::string&);
        
void read_qual_mapper::map(const std::string& key, const std::string& value, kv_sink& context)
{
	read_id id;
	unaligned_reads reads;
	try
	{
		msgpack_deserialize(id, key);
		msgpack_deserialize(reads, value);
	}
	catch(const deserialization_error& e)
	{
		throw io_exception(printstring("in read_qual_mapper: %s", e.message().c_str())); 
	}
	
	std::string one = msgpack_serialize(1UL);
	std::string buf(3,0);

	for(size_t i = 0; i < reads.size(); i++)
	{
		buf.resize(3);
		for(size_t j = 0; j < reads[i].sequence.size(); j++)
		{
			buf[0] = reads[i].sequence[j];
			buf[1] = reads[i].quality[j] - 33;
			buf[2] = j + 1;
			context.write(buf, one);
		}
		buf.resize(2);
		buf[0] = 'E';
		buf[1] = reads[i].sequence.size() + 1;
		context.write(buf, one);
	}
}


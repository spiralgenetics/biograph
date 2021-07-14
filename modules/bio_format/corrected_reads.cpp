#include "modules/bio_format/corrected_reads.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/io/log.h"
#include "modules/io/registry.h"

#include <boost/format.hpp>


REGISTER_3(exporter, corrected_reads, writable&, bool, const std::string&);

void corrected_reads_exporter::write(const std::string& key, const std::string& value)
{
	std::string name;
	corrected_reads the_reads;
	msgpack_deserialize(name, key);
	msgpack_deserialize(the_reads, value);
	for(size_t i = 0; i < the_reads.size(); i++)
	{
		const corrected_read& the_read = the_reads[i];
		std::string a_line =
			boost::str(boost::format("%1%\t%2%\t%3%\t%4%\t%5%\n")
			//boost::str(boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\n")
			% name 
			% the_read.sequence.as_string() 
			% the_read.corrected.as_string() 
			% the_read.quality
			% the_read.trace_me
			)
		;
		m_sink.write(a_line.c_str(), a_line.size());
		//m_sink.print("@%s\n%s\n+\n%s\n", name.c_str(), the_read.corrected.as_string().c_str(), the_read.quality.c_str());
	}
}

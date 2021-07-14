
#include <boost/foreach.hpp>
#include <map>
#include <string.h>

#include "modules/bio_base/aligned_read.h"
#include "modules/bio_base/karyotype_compat.h"
#include "modules/bio_base/seq_position.h"
#include "modules/mapred/ex_im_porter_data.h"
#include "modules/bio_format/sam_type.h"
#include "modules/io/log.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/registry.h"

REGISTER_3(importer, sam, readable&, bool, const std::string&);
REGISTER_3(exporter, sam, writable&, bool, const std::string&);


sam_exporter::sam_exporter(writable & byte_sink, const std::string& ref_name_name, bool use_supercontig_coords,
	const std::string& start_key, const std::string& end_key)
	: exporter(byte_sink), m_reference(NULL), m_use_supercontig_coords(use_supercontig_coords)
{
	set_ref_name(m_reference, ref_name_name);
	m_are_reads_sorted = (! start_key.empty()) && (! end_key.empty());
	if (m_are_reads_sorted)
	{
		msgpack_deserialize(m_start_read_pos, start_key);
		msgpack_deserialize(m_end_read_pos, end_key);
	}
}

sam_exporter::sam_exporter(writable & byte_sink,
	bool /*serialize_argument*/,
	const std::string& serialized_exporter_data)
	: exporter(byte_sink), m_reference(NULL)
{
	if (! serialized_exporter_data.empty())
	{
		ex_im_porter_data the_exporter_data;
		msgpack_deserialize(the_exporter_data, serialized_exporter_data);
		set_ref_name(m_reference, the_exporter_data.ref_name);
		m_use_supercontig_coords = false;

		m_are_reads_sorted = (! the_exporter_data.start_key.empty()) && (! the_exporter_data.end_key.empty());
		if (m_are_reads_sorted)
		{
			msgpack_deserialize(m_start_read_pos, the_exporter_data.start_key);
			msgpack_deserialize(m_end_read_pos, the_exporter_data.end_key);
		}
	}
}

void sam_exporter::write_header()
{
	if (m_reference != NULL)
	{
		m_sink.print("@HD\tVN:1.3\tSO:coordinate\n");

		reference_assembly ref = m_reference->get_assembly();

		for(int i = 0; i < (int) ref.scaffold_order.size(); i++)
		{
			scaffold chrom = ref.get_scaffold(ref.scaffold_order[i]);
			// if ((! m_are_reads_sorted) || (the_karyotype.does_range_contain_chromosome(m_start_read_pos, m_end_read_pos, name)))
			// {
				m_sink.print("@SQ\tSN:%s\tLN:%ld\n", chrom.name.c_str(), chrom.len);
			// }
		}

		// HACK:  Print a fake read group to make GATK happy.  Remove when
		// metadata handling is implemented.
		m_sink.print("@RG\tID:Spiral\tSM:Sample1\tLB:Library\tPL:Illumina\n");
	}
}

void sam_exporter::write(const std::string& /*key*/, const std::string& value)
{
	aligned_read an_aligned_read;
	msgpack_deserialize(an_aligned_read, value);

	m_sink.print("%s\n", print_sam(m_reference->get_assembly(), an_aligned_read, m_use_supercontig_coords).c_str());
}

sam_importer::sam_importer(readable& source, bool /*unused*/, const std::string& serialized_importer_data)
	: m_source(source)
{
	ex_im_porter_data the_importer_data;
	msgpack_deserialize(the_importer_data, serialized_importer_data);
	set_ref_name(m_reference, the_importer_data.ref_name);
}

sam_importer::sam_importer(readable& source, const std::string& ref_name_name) : m_source(source)
{
	set_ref_name(m_reference, ref_name_name);
}

void sam_importer::import(kv_sink& sink, simple_metadata& meta)
{
	SPLOG("Importing SAM");
	std::string	sam_line;
	while(m_source.readline(sam_line, 4096))
	{
		if (sam_line[0] == '@')
		{
			continue;
		}

		aligned_read the_read;
		if (parse_sam(m_reference->get_assembly(), the_read, sam_line))
		{
			if (!(the_read.flags & 4))
			{
				sink.write_msgpack(the_read.ref_pos, the_read);
			}
		}
		else
		{
			SPLOG("Unable to import sam: %s", sam_line.c_str());
			break;
		}
	}
	SPLOG("Done importing SAM");
}

void set_ref_name(boost::scoped_ptr<reference>& a_reference, const std::string& ref_name_name)
{
	if (! ref_name_name.empty())
		a_reference.reset(new reference(ref_name_name));
}


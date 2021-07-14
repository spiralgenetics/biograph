
#include "modules/bio_format/assembly.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_format/struct_var.h"
#include "modules/io/log.h"
#include "modules/io/registry.h"

#include <boost/format.hpp>

#include <string.h>
#include <map>

REGISTER_3(exporter, assembly, writable&, bool, const std::string&);

assembly_exporter::assembly_exporter(
	writable& byte_sink,
	const std::string& ref_name)
	: exporter(byte_sink)
	, m_last_var_id(-1)
	, m_reference(ref_name)
	, m_reference_assembly(m_reference.get_assembly())
{}

assembly_exporter::assembly_exporter(
	writable& byte_sink,
	bool /*unused*/,
	const std::string& ref_name)
	: exporter(byte_sink)
	, m_last_var_id(-1)
	, m_reference(ref_name)
	, m_reference_assembly(m_reference.get_assembly())
{}

void assembly_exporter::write(const std::string & key, const std::string & value)
{
	struct_var sv;
	msgpack_deserialize(sv, value);

	// TODO: don't filter by SV
	// TODO: add vartype

	// format is as follows:
	// Assembly ID
	// Assembled Sequence
	// Reference Location 1 ([+-]SeqName:SeqOffset)
	// Reference Location 2 ([+-]SeqName:SeqOffset)
	// Breakend Offset 1
	// Breakend Offset 2
	// Edit Distance
	// Depth
	const char* FORMAT = "%d\t%s\t%c%s:%d\t%c%s:%d\t%d\t%d\t%d\t%d\n";

	std::string line = boost::str(boost::format(FORMAT)
		% sv.var_id
		% sv.assembled.as_string()
		% (sv.rev_start ? '-' : '+')
		% m_reference_assembly.scaffold_order[sv.ref_start.scaffold_id]
		% sv.ref_start.position
		% (sv.rev_end ? '-' : '+')
		% m_reference_assembly.scaffold_order[sv.ref_end.scaffold_id]
		% sv.ref_end.position
		% sv.var_start
		% sv.var_end
		% sv_compute_edit_distance(sv, m_reference)
		% sv.depth
	);
	m_sink.write(line.c_str(), line.size());

	m_last_var_id = sv.var_id;
}

#pragma once

#include "modules/bio_base/reference.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_base/align_astar.h"
#include "modules/io/config.h"
#include "modules/bio_format/exporter.h"

int sv_compute_edit_distance(const struct_var& sv, const reference& ref);

class struct_var_exporter : public exporter
{
public:
	struct_var_exporter(writable & byte_sink, const std::string& ref_name)
		: exporter(byte_sink)
		, m_reference(ref_name)
		, m_reference_assembly(m_reference.get_assembly())
	{}

	struct_var_exporter(writable & byte_sink,
		bool /*unused*/,
		const std::string& ref_name)
		: exporter(byte_sink)
		, m_reference(ref_name)
		, m_reference_assembly(m_reference.get_assembly())
	{}

	void write(const std::string& key, const std::string& value) override;

private:
	reference m_reference;
	const reference_assembly& m_reference_assembly;
};

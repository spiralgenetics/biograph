#pragma once

#include "modules/bio_format/exporter.h"
#include "modules/bio_base/reference.h"

class assembly_exporter : public exporter
{
public:
	assembly_exporter(
		writable& byte_sink,
		const std::string& ref_name);

	assembly_exporter(
		writable& byte_sink,
		bool /*unused*/,
		const std::string& ref_name);

public:
	void write(const std::string& key, const std::string& value) override;

private:
	int m_last_var_id;
	reference m_reference;
	const reference_assembly& m_reference_assembly;
};

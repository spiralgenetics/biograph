#pragma once

#include "modules/bio_base/reference.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_base/align_astar.h"
#include "modules/io/registry.h"
#include "modules/bio_format/exporter.h"
#include "modules/mapred/path.h"
#include "modules/io/config.h"

struct sqlite3;
struct sqlite3_stmt;

class svdb_exporter : public exporter
{
public:
	svdb_exporter(writable& byte_sink, const std::string& ref_name);

	svdb_exporter(writable& byte_sink,
		bool /*unused*/,
		const std::string& ref_name)
		: exporter(byte_sink)
		, m_reference(ref_name)
		, m_reference_assembly(m_reference.get_assembly())
		, m_closed(false)
	{}

	virtual void write(const std::string& key, const std::string& value);
	virtual void close();

private:
	reference m_reference;
	const assembly& m_reference_assembly;
	sqlite3 *m_db;
	sqlite3_stmt *m_insert_prepared_statement;
	bool m_closed;

};

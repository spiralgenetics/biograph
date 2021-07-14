#pragma once

#include <boost/scoped_ptr.hpp>

#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/io/keyvalue.h"
#include "modules/bio_format/importer.h"
#include "modules/bio_format/exporter.h"

class serialized_pipe_params;
class sam_exporter : public exporter
{
public:
	sam_exporter(writable & byte_sink, const std::string& ref_name, bool use_supercontig_coords = false,
	// If the start and end keys are omitted, the reads are assumed unsorted and the entire species karyotype
	// is written to the header.  The keys should be magpack serialized.
	const std::string& start_key = std::string(), const std::string& end_key = std::string());
	sam_exporter(writable & byte_sink,
		bool /*serialize_argument*/, // Used to disambiguate the constructors.  The actual value is unused.
		const std::string& serialized_pipe_params);

	void write_header() override;
	void write(const std::string & key, const std::string & value) override;


private:
	boost::scoped_ptr<reference> m_reference;
	seq_position m_start_read_pos;
	seq_position m_end_read_pos;
	bool m_use_supercontig_coords;
	bool m_are_reads_sorted;
};

class sam_importer : public importer
{
public:
	sam_importer(readable& source, bool /*unused*/, const std::string& serialized_pipe_params);
	sam_importer(readable& source, const std::string& ref_name);

  void import(kv_sink& sink, simple_metadata& meta) override;


private:
	readable& m_source;
	boost::scoped_ptr<reference> m_reference;
};

void set_ref_name(boost::scoped_ptr<reference>& a_reference, const std::string& ref_name);


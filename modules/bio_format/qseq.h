#pragma once

#include "modules/bio_format/importer.h"

class qseq_importer : public importer
{
public:
	qseq_importer(readable& source) 
		: m_source(source) 
	{}
	qseq_importer(readable& source, bool /*unused*/, const std::string& /*unused_pipe_params*/)
		: m_source(source) 
	{}

	void import(kv_sink& sink, simple_metadata& meta) override;


private:
	readable& m_source;
};

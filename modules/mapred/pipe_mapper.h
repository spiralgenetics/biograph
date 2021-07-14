#pragma once

#include <deque>

#include "modules/io/keyvalue.h"
#include "modules/bio_format/exporter.h"
#include "modules/mapred/map_pipe_task.h"

class unix_pipeline;
class pipe_mapper_buffer : public readable, public writable
{
	friend class pipe_mapper;

public:
	pipe_mapper_buffer(kv_source& kv_data_source, map_pipe_task* const the_map_pipe_task)
		: m_exporter(nullptr)
		, m_unix_pipeline(nullptr)
		, m_source(kv_data_source)
		, m_closed(false)
		, m_map_pipe_task(the_map_pipe_task) {}

	size_t read(char* buf, size_t len) override;
	void write(const char* buf, size_t len) override;

private:
	// This is a godawful hack.  TODO: fix this abomination.
	void set_exporter(exporter& an_exporter) {
		m_exporter = &an_exporter;
		m_exporter->write_header();
	}
	exporter* m_exporter;

	void set_unix_pipeline(unix_pipeline& a_unix_pipeline) {
		m_unix_pipeline = &a_unix_pipeline;
	}
	unix_pipeline* m_unix_pipeline;

	kv_source& m_source;
	std::deque<char> m_buffer;
	bool m_closed;
	map_pipe_task* const m_map_pipe_task;
};


// To use a pipe_mapper you need to supply an importer suitable for the kv_data
// you want to write along with an exporter for the same type of kv data. The
// importer and exporter will require a readable and writable respectively.
// The readable for the importer is the buffer class and the writable for the
// exporter is a unix_pipeline object that will need to be created and passed
// in. The unix_pipeline target writable is the buffer. Once this is done,
// pass the importer, exporter and buffer into a new pipe_mapper. Then when
// you're ready to map some data, call operator() on your pipe_mapper object

class pipe_mapper
{
public:
	pipe_mapper(
		pipe_mapper_buffer& a_buffer,
		exporter& an_exporter,
		importer& an_importer,
		unix_pipeline& a_unix_pipeline
	) : m_buffer(a_buffer)
		, m_importer(an_importer)
		{
			m_buffer.set_exporter(an_exporter);
			m_buffer.set_unix_pipeline(a_unix_pipeline);
		}

	void operator() (kv_sink& kv_data_sink) {
		m_importer.import(kv_data_sink);
	}

private:
	pipe_mapper_buffer& m_buffer;
	importer& m_importer;
};

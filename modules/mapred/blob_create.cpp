
#include "modules/mapred/blob_create.h"
#include "modules/io/uuid.h"

blob_create::blob_create(manifest& out, const path& prefix, const std::string& job_id, size_t max_in_flight, size_t chunk_size)
	: m_out(out)
	, m_prefix(prefix)
	, m_job_id(job_id)
	, m_max_in_flight(max_in_flight)
	, m_chunk_size(chunk_size)
{
	m_in_flight.emplace_back();
}

void blob_create::write(const char* buf, size_t len) {
	while(len) {
		if (m_in_flight.back().buffer.size() == m_chunk_size) {
			finish_chunk();
			m_in_flight.emplace_back();
		}
		pending_chunk& chunk = m_in_flight.back();
		size_t to_write = std::min(len, m_chunk_size - chunk.buffer.size());
		chunk.buffer.write(buf, to_write);
		buf += to_write;
		len -= to_write;
	}
}

void blob_create::close() {
	if (m_in_flight.back().buffer.size() == 0) {
		m_in_flight.pop_back();
	} else {
		finish_chunk();
	}
	while(!m_in_flight.empty()) {
		m_in_flight.front().wait->wait();
		m_in_flight.pop_front();
	}
}

void blob_create::finish_chunk() {
	pending_chunk& chunk = m_in_flight.back();
	chunk.write_path = m_prefix.append(printstring("%03ld_%s",
		random() % 1000,
		make_uuid().c_str()));

	path_write_options options;
	if (!m_job_id.empty()) {
		options.tags["Job"] = m_job_id;
	}
	chunk.wait = chunk.write_path.write_inverted(chunk.buffer, chunk.buffer.size(), options);
	file_info fi(chunk.write_path, chunk.buffer.size(), 0);
	m_out.add(fi, 0);
	if (m_in_flight.size() == m_max_in_flight) {
		m_in_flight.front().wait->wait();
		m_in_flight.pop_front();
	}
}

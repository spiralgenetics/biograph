#include <algorithm>
#include <cassert>
#include <string>
#include <tr1/array>
#include <unistd.h>
#include <vector>

#include "base/base.h"
#include "modules/io/log.h"
#include "modules/bio_format/exporter.h"
#include "modules/bio_format/importer.h"
#include "modules/mapred/pipe_mapper.h"
#include "modules/mapred/pipe_params.h"
#include "modules/mapred/unix_pipeline.h"

size_t pipe_mapper_buffer::read(char* buf, size_t len)
{
	// SPLOG("pipe_mapper_buffer::read> reading, m_closed = %s, m_buffer.size() = %lu, len = %lu", m_closed ? "true" : "false", m_buffer.size(), len);

	while (len > m_buffer.size() && m_unix_pipeline->is_child_alive()) {
		if (m_closed) {
			usleep(100000);
			continue;
		}

		std::string key;
		std::string value;
		if (m_source.read(key, value)) {
			CHECK(m_exporter);
			m_exporter->write(key, value);
			if (m_map_pipe_task) {
				m_map_pipe_task->processed_a_record();
			}
		}
		else {
			m_exporter->write_footer();
			m_exporter->close();
			m_closed = true;
		}
	}

	if (m_closed && len > m_buffer.size()) {
		len = m_buffer.size();
	}
	std::copy(m_buffer.rbegin(), m_buffer.rbegin() + len, buf);
	m_buffer.erase((m_buffer.rbegin() + len).base(), m_buffer.rbegin().base());

	return len;
}


void pipe_mapper_buffer::write(const char* buf, size_t len)
{
	std::copy(buf, buf + len, front_inserter(m_buffer));
}


#pragma once

#include <list> 

#include "modules/io/mem_io.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/path.h"

class blob_create : public writable {
public:
	blob_create(manifest& out, const path& prefix, const std::string& job_id /* otherwise known as "user" or "project" */, size_t max_in_flight = 8, size_t chunk_size = 64*1024*1024);
	void write(const char* buf, size_t len) override;
	void close() override;
private:
	void finish_chunk();
	struct pending_chunk {
      mem_io buffer{"",track_alloc("blob_create:pending_chunk")};
		path write_path;
		std::unique_ptr<waiter> wait;
	};
	manifest& m_out;
	path m_prefix;
	std::string m_job_id;
	size_t m_max_in_flight;
	size_t m_chunk_size;
	std::list<pending_chunk> m_in_flight;
};


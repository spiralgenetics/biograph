#pragma once

#include "modules/io/mmap_buffer.h"
#include "modules/mapred/manifest.h"
#include "modules/io/progress.h"

class resource_manager
{
public:
	resource_manager();

	// this constructor is used for unit testing
	explicit resource_manager(bool direct);

	void create_resource(mmap_buffer& out, size_t size);

	void write_resource(
		manifest& out,
		mmap_buffer& in,
		const path& root,
		const std::string& prefix,
		const progress_handler_t& progress = null_progress_handler
	);

	void read_resource(
		mmap_buffer& out,
		const manifest& in,
		const progress_handler_t& progress = null_progress_handler
	);

private:
	size_t free_space();
	void make_space(size_t size);
	bool m_direct = false;
};

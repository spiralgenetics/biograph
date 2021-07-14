#pragma once

#include "modules/io/progress_tracker_types.h"

class progress_tracker
{
	public:
		progress_tracker(progress_t& update);

		void update(size_t in, size_t out);
		void final_update(size_t in, size_t out);

	private:
		progress_t	m_update;
		size_t		m_progress_chunk_size;
		size_t		m_progress_chunks_completed;
};

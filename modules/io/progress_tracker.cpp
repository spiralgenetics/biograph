#include "modules/io/progress_tracker.h"

	progress_tracker::progress_tracker(progress_t& update)
:m_update(update),
	m_progress_chunk_size(1),
	m_progress_chunks_completed(0)
{ }

void progress_tracker::update(size_t in, size_t out)
{
	if ( m_progress_chunks_completed != (in / m_progress_chunk_size) )
	{
		m_progress_chunk_size = m_update(in, out);
		// re-compute current progress, in case m_update returned
		// a new modulo.
		m_progress_chunks_completed = in / m_progress_chunk_size;
	}
}

void progress_tracker::final_update(size_t in, size_t out)
{
	(void) m_update(in, out);
}

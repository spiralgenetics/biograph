#pragma once
#include <mutex>
#include <condition_variable>

#include "modules/io/io.h"
#include "modules/io/log.h"

// This class acts as a barrier to synchronize multiple threads so that they
// all wait until everyone is done and then proceed.  To use, initialize the
// class with the number of threads, then each thread calls the wait method.
// When the first through N-1th threads call wait, they will block until the
// Nth thread calls wait.   At that point all threads will wake and proceed.
class thread_barrier
{
public:
	explicit thread_barrier(int count)
		: m_thread_count(count)
		, m_original_thread_count(count)
		, m_generation(0)
	{}

	void wait()
	{
		std::unique_lock<std::mutex> lock{ m_mutex };
		SPLOG("thread_barrier::wait()> count = %d, generation = %d", m_thread_count, m_generation);
		int local_generation = m_generation;
		
		if (--m_thread_count == 0) {
			SPLOG("Last thread!");
			m_generation++;
			m_thread_count = m_original_thread_count;
			m_cv.notify_all();
			return;
		}

		m_cv.wait(
			lock
			, [this, local_generation] { return local_generation != m_generation; }
		);
	}
	
	void reset(int thread_count)
	{
		if (m_thread_count) {
			throw io_exception("thread_barrier::reset> Reset called when threads are waiting: "
				+ std::to_string(m_thread_count));
		}
		
		m_thread_count = thread_count;
	}
	
private:
	std::mutex m_mutex;
	std::condition_variable m_cv;
	int m_thread_count;
	int m_original_thread_count;
	int m_generation;
};

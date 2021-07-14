#pragma once

#include <mutex>
#include <condition_variable>

class semaphore
{
public:
	semaphore()
		: m_count(0)
	{}

	explicit semaphore(size_t count)
		: m_count(count)
	{}

	void notify()
	{
		std::lock_guard<std::mutex> lock(m_mutex);
		++m_count;
		m_condition.notify_one();
	}

	void wait()
	{
		std::unique_lock<std::mutex> lock(m_mutex);
		while (!m_count) {
			m_condition.wait(lock);
		}
		--m_count;
	}

private:
	std::mutex m_mutex;
	std::condition_variable m_condition;
	size_t m_count;
};


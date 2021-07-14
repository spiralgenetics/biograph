#pragma once

#include <thread>
#include <functional>
#include <condition_variable>

typedef std::function<void(std::cv_status)> pulse_f;

class pulse
{
public:
	pulse(std::chrono::seconds period, pulse_f fn)
		: m_period(period)
		, m_fn(fn)
		, m_thread(&pulse::run, this)
	{}

	~pulse()
	{
		stop();
	}

	void stop()
	{
		if (m_terminate) {
			return;
		}

		{
			std::lock_guard<std::mutex> lock(m_mutex);
			m_terminate = true;
		}
		m_cond.notify_all();
		m_thread.join();
	}

private:
	void run()
	{
		std::unique_lock<std::mutex> lock(m_mutex);
		while (!m_terminate) {
			// timed_wait will block until either the condition is set (by a signal handler) 
			// or period has elapsed.
			auto status = m_cond.wait_for(lock, m_period);
			try {
				m_fn(status);
			}
			catch (const std::exception& ex) {
				SPLOG("pulse::run> uncaught exception: %s", ex.what());
			}
		}
	}

private:
	volatile bool              m_terminate = false;
	std::chrono::seconds       m_period;
	pulse_f                    m_fn;
	std::mutex                 m_mutex;
	std::condition_variable    m_cond;
	std::thread                m_thread;
};

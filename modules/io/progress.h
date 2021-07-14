#pragma once

#include "base/base.h"
#include "modules/io/log.h"
#include <functional>
#include <future>

typedef std::function<void(double)> progress_handler_t;

inline void null_progress_handler(double x) {}
inline void null_watchdog() {}

class subprogress
{
public:
	inline subprogress(const progress_handler_t& outer, double start, double end)
		: m_outer(outer)
		, m_start(start)
		, m_end(end) 
	{}

	inline void operator()(double progress) 
	{ 
		m_outer(m_start + progress * (m_end - m_start)); 
	}

private:
	progress_handler_t m_outer;
	double m_start;
	double m_end;
};

// A version of subprogress that divides the progress into equal parts.
class equal_subprogress {
 public:
  equal_subprogress(const progress_handler_t& outer, size_t num_parts)
      : m_outer(outer), m_num_parts(num_parts) {}

  progress_handler_t operator[](size_t part_num) const {
    CHECK_LT(part_num, m_num_parts);
    return subprogress(m_outer, part_num * 1. / m_num_parts,
                       (part_num + 1) * 1. / m_num_parts);
  }

 private:
  progress_handler_t m_outer;
  size_t m_num_parts;
};

class noisy_progress_handler
{
public:
	inline void operator()(double progress) 
	{
		time_t now = time(0);
		if (m_last_update + 2 <= now) { 
			m_last_update = now;
			SPLOG("Progress: %f", progress);
		}
	}

private:
	time_t m_last_update = 0;
};

template <typename T>
T future_watchdog(
	std::future<T>& future, 
	const progress_handler_t& on_progress, 
	double progress)
{
	while (true) {
		auto status = future.wait_for(std::chrono::seconds(1));
		if (status == std::future_status::ready) {
			break;
		}
		on_progress(progress);
	}
	return future.get();
}

template <typename F>
void lambda_watchdog(const progress_handler_t& on_progress, double progress, const F& func) {
	std::future<void> f = std::async(std::launch::async, func);
	future_watchdog(f, on_progress, progress);
}


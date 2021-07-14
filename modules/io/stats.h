#pragma once

#include <vector>
#include <limits>

template <typename T>
struct simple_stats
{
	simple_stats()
		: avg(0)
		, min(std::numeric_limits<T>::max())
		, max(0)
	{}

	void add_sample(const T& value)
	{
		samples.push_back(value);
	}

	void analyze()
	{
		double sum = 0.0;
		for (const auto& sample : samples) {
			min = std::min(min, sample);
			max = std::max(max, sample);
			sum += sample;
		}
		avg = sum / samples.size();
	}

	std::vector<T> samples;
	double avg;
	T min;
	T max;
};

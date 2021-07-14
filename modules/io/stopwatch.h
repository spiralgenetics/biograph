#pragma once

#include <chrono>

template <typename Function>
std::chrono::milliseconds stopwatch(const Function& fn) {
	auto start = std::chrono::high_resolution_clock::now();
	fn();
	auto stop = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
}

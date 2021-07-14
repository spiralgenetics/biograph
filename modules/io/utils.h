#pragma once

#include <cstring>
#include <ctime>
#include <string>

#include <boost/filesystem.hpp>

#include "base/base.h"

namespace fs = boost::filesystem;

// Helper to 'printf' to strings
std::string printstring(const char* format, ...) __attribute__((format(printf,1,2)));

template <class Function>
void map_lines(const std::string& str, const Function& fn)
{
	size_t pos = 0;
	while (pos < str.size()) {
		auto next = str.find('\n', pos);
		if (next == std::string::npos) {
			fn(str.substr(pos));
			break;
		}
		fn(str.substr(pos, next - pos));
		pos = next + 1;
	}
}

std::string time_to_RFC3339(const std::time_t time);

// Print a progress bar to the screen (0..1)
void print_progress(float progress, int width = 70);

// Enable/disable terminal echo
void setecho(bool enable);

// Number of columns in the current terminal
unsigned int get_terminal_width();

// Set a strict memory limit (NOTE: also counts against mmaps!)
void set_mem_limit(uint64_t max_mem);

// Get the current memory limit
uint64_t get_mem_limit();

// Total system memory
uint64_t get_system_mem();

// kernel info
std::string get_uname();

// node name
std::string get_nodename();

// Linux release
std::string get_os_release();

//Strongly typed memcpys to avoid repeating oneself
template<typename T> void* typed_memcpy(void* dest, const T* source)
{
	return std::memcpy(dest, source, sizeof(T));
}

// Also works for std::strings or anything with contiguous storage, size method
// and value_type typedef.
template<typename V> void* vector_memcpy(void* dest, const V& a_vector)
{
	return std::memcpy(dest, a_vector.data(), a_vector.size() * sizeof(typename V::value_type));
}

// Run a command through shell, return the output as a std::string
std::string easy_exec(const char* cmd);
std::string easy_exec(std::string cmd);

// Gets the minimum value from a container.
template <typename T>
auto container_min(const T& container) -> typename T::value_type {
  CHECK(!container.empty());
  auto it = container.begin();
  CHECK(it != container.end());
  typename T::value_type result = *it;
  ++it;
  while (it != container.end()) {
    if (*it < result) {
      result = *it;
    }
    ++it;
  }
  return result;
}

std::string expand_home(const std::string& path);

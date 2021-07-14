#pragma once

// C++17-style string_view class.
//
// When we move to C++17, we should be able to replace this with
// std::string_view.

#include <iostream>

class string_view {
 public:
  string_view() = default;
  string_view(const std::string& str)
      : m_data(str.data()), m_size(str.size()) {}
  string_view(const char* data, size_t size) : m_data(data), m_size(size) {}
  string_view(const char* data) : m_data(data), m_size(strlen(data)) {}

  const char* begin() const { return m_data; }
  const char* end() const { return m_data + m_size; }

  operator std::string() const { return std::string(m_data, m_size); }

  const char* data() const { return m_data; }
  size_t size() const { return m_size; }

 private:
  const char* m_data = nullptr;
  size_t m_size = 0;
};

inline std::ostream& operator<<(std::ostream& os, const string_view& sv) {
  return os << std::string(sv);
}

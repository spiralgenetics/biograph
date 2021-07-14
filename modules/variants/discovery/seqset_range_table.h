#pragma once

#include "modules/bio_base/seqset.h"

#include <boost/range.hpp>

namespace variants {
namespace discovery {

struct seqset_range_comparer {
  bool operator()(const seqset_range& lhs, const seqset_range& rhs) const {
    if (lhs.begin() != rhs.begin()) {
      return lhs.begin() < rhs.begin();
    }
    // Less specific before more specific.
    return lhs.size() < rhs.size();
  }
};

template <typename V>
class seqset_range_table {
  using data_t = std::map<seqset_range, V, seqset_range_comparer>;

 public:
  seqset_range_table() = default;
  seqset_range_table(const seqset_range_table&) = delete;
  seqset_range_table& operator=(const seqset_range_table&) = delete;

  using iterator = typename data_t::iterator;
  using const_iterator = typename data_t::const_iterator;

  iterator begin() { return m_data.begin(); }
  iterator end() { return m_data.end(); }
  const_iterator begin() const { return m_data.begin(); }
  const_iterator end() const { return m_data.end(); }

  using iterator_range = boost::iterator_range<iterator>;
  iterator_range entries_starting_with(const seqset_range& r);
  using const_iterator_range = boost::iterator_range<const_iterator>;
  const_iterator_range entries_starting_with(const seqset_range& r) const;

  V& operator[](const seqset_range& r) { return m_data[r]; }
  V& at(const seqset_range& r) { return m_data[r]; }
  const V& at(const seqset_range& r) const { return m_data.at(r); }

 private:
  data_t m_data;
};

template <typename V>
typename seqset_range_table<V>::iterator_range seqset_range_table<V>::entries_starting_with(
    const seqset_range& r) {
  iterator begin = m_data.lower_bound(r);
  iterator end = begin;
  while (end != m_data.end() && end->first.end() <= r.end()) {
    ++end;
  }
  return iterator_range(begin, end);
}
template <typename V>
typename seqset_range_table<V>::const_iterator_range seqset_range_table<V>::entries_starting_with(
    const seqset_range& r) const {
  const_iterator begin = m_data.lower_bound(r);
  const_iterator end = begin;
  while (end != m_data.end() && end->first.end() <= r.end()) {
    ++end;
  }
  return const_iterator_range(begin, end);
}

}  // namespace discovery
}  // namespace variants

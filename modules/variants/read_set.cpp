#include "modules/variants/read_set.h"

namespace variants {

void read_coverage_set::insert(int offset, uint32_t read_id, int read_len) {
  len_and_offset loff;
  loff.len = read_len;
  loff.offset = offset;

  auto& reads = m_impl[loff];
  reads.insert(read_id);
}

void read_coverage_set::insert(const read_coverage_read_t& new_cov) {
  len_and_offset loff;
  loff.len = new_cov.read_len;
  loff.offset = new_cov.offset;

  auto& reads = m_impl[loff];

  reads.insert(new_cov.read_ids);
}

void read_coverage_set::insert(int offset, const read_id_set& read_ids, int read_len) {
  len_and_offset loff;
  loff.len = read_len;
  loff.offset = offset;

  auto& reads = m_impl[loff];

  reads.insert(read_ids);
}

constexpr size_t read_id_set::k_mask_bits;
constexpr size_t big_read_id_set::k_mask_bits;
constexpr size_t read_id_set::k_num_small_elem;

read_id_set::read_id_set(const big_read_id_set& orig) {
  m_impl.reserve(orig.m_impl.size());
  for (const auto& e : orig.m_impl) {
    elem new_elem;
    new_elem.chunk_id = e.first;
    new_elem.read_id_bits = e.second;
    m_impl.push_back(new_elem);
  }
}

size_t read_id_set::size() const {
  size_t tot = 0;
  for (const auto& e : m_impl) {
    tot += __builtin_popcountl(e.read_id_bits);
  }
  return tot;
}

std::ostream& operator<<(std::ostream& os, const read_id_set& ids) {
  bool first = true;
  os << "ReadIds(";
  for (uint32_t read_id : ids) {
    if (first) {
      first = false;
    } else {
      os << ",";
    }
    os << read_id;
  }
  return os << ")";
}

void read_id_set::insert(const read_id_set& old_ids) {
  if (m_impl.empty()) {
    m_impl = old_ids.m_impl;
    return;
  }

  for (const auto& old_elem : old_ids.m_impl) {
    auto it = std::lower_bound(m_impl.begin(), m_impl.end(), old_elem.chunk_id);

    if (it != m_impl.end() && it->chunk_id == old_elem.chunk_id) {
      it->read_id_bits |= old_elem.read_id_bits;
    } else {
      if (it != m_impl.end()) {
        DCHECK_GT(it->chunk_id, old_elem.chunk_id);
      }
      m_impl.insert(it, old_elem);
    }
  }
}

void read_id_set::insert(uint32_t read_id) {
  uint32_t chunk_id = read_id / k_mask_bits;
  unsigned offset = read_id & (k_mask_bits - 1);
  auto it = std::lower_bound(m_impl.begin(), m_impl.end(), chunk_id);
  if (it != m_impl.end() && it->chunk_id == chunk_id) {
    it->read_id_bits |= 1ULL << offset;
    return;
  }

  elem rs;
  rs.chunk_id = chunk_id;
  rs.read_id_bits |= 1UL << offset;
  m_impl.insert(it, std::move(rs));
}

void read_id_set::erase(uint32_t read_id) {
  uint32_t chunk_id = read_id / k_mask_bits;
  unsigned offset = read_id & (k_mask_bits - 1);

  auto it = std::lower_bound(m_impl.begin(), m_impl.end(), chunk_id);
  if (it != m_impl.end() && it->chunk_id == chunk_id) {
    it->read_id_bits &= ~(1ULL << offset);

    if (it->read_id_bits == 0) {
      m_impl.erase(it);
    }
  }
}

void read_id_set::clear() { m_impl.clear(); }

bool read_id_set::contains(uint32_t read_id) const {
  uint32_t chunk_id = read_id / k_mask_bits;
  unsigned offset = read_id & (k_mask_bits - 1);

  auto it = std::lower_bound(m_impl.begin(), m_impl.end(), chunk_id);
  if (it != m_impl.end() && it->chunk_id == chunk_id) {
    return it->read_id_bits & (1ULL << offset);
  }
  return false;
}

read_id_set read_id_set::intersection(const read_id_set& rhs) const {
  read_id_set result;
  result.m_impl.reserve(std::min(m_impl.size(), rhs.m_impl.size()));

  auto it1 = m_impl.begin();
  auto it2 = rhs.m_impl.begin();
  while (it1 != m_impl.end() && it2 != rhs.m_impl.end()) {
    if (it1->chunk_id < it2->chunk_id) {
      ++it1;
      continue;
    }
    if (it2->chunk_id < it1->chunk_id) {
      ++it2;
      continue;
    }

    CHECK_EQ(it1->chunk_id, it2->chunk_id);
    elem intersection = *it1;
    intersection.read_id_bits &= it2->read_id_bits;

    if (!intersection.read_id_bits) {
      ++it1;
      ++it2;
      continue;
    }

    result.m_impl.emplace_back(std::move(intersection));

    ++it1;
    ++it2;
  }
  return result;
}

std::vector<uint32_t> read_id_set::to_vector() const {
  std::vector<uint32_t> result(begin(), end());
  return result;
}

read_id_set read_id_set::operator|(const read_id_set& rhs) const {
  read_id_set result;

  result.m_impl.reserve(m_impl.size() + rhs.m_impl.size());

  auto it = m_impl.begin();
  auto rhs_it = rhs.m_impl.begin();

  while (it != m_impl.end() && rhs_it != rhs.m_impl.end()) {
    if (it->chunk_id < rhs_it->chunk_id) {
      result.m_impl.emplace_back(*it);
      ++it;
      continue;
    }
    if (rhs_it->chunk_id < it->chunk_id) {
      result.m_impl.emplace_back(*rhs_it);
      ++rhs_it;
      continue;
    }

    CHECK_EQ(it->chunk_id, rhs_it->chunk_id);
    result.m_impl.emplace_back(*it);
    result.m_impl.back().read_id_bits |= rhs_it->read_id_bits;
    ++it;
    ++rhs_it;
  }

  std::copy(it, m_impl.end(), std::back_inserter(result.m_impl));
  std::copy(rhs_it, rhs.m_impl.end(), std::back_inserter(result.m_impl));

  return result;
}

read_id_set read_id_set::operator-(const read_id_set& rhs) const {
  read_id_set result;

  result.m_impl.reserve(m_impl.size());

  auto it = m_impl.begin();
  auto rhs_it = rhs.m_impl.begin();

  while (it != m_impl.end() && rhs_it != rhs.m_impl.end()) {
    if (it->chunk_id < rhs_it->chunk_id) {
      result.m_impl.emplace_back(*it);
      ++it;
      continue;
    }
    if (rhs_it->chunk_id < it->chunk_id) {
      ++rhs_it;
      continue;
    }

    CHECK_EQ(it->chunk_id, rhs_it->chunk_id);
    elem new_elem;
    new_elem.chunk_id = it->chunk_id;
    new_elem.read_id_bits = it->read_id_bits & ~rhs_it->read_id_bits;
    if (new_elem.read_id_bits) {
      result.m_impl.emplace_back(new_elem);
    }
    ++it;
    ++rhs_it;
  }

  std::copy(it, m_impl.end(), std::back_inserter(result.m_impl));

  return result;
}

read_id_set read_id_set::operator&(const read_id_set& rhs) const {
  read_id_set result;

  result.m_impl.reserve(std::min(m_impl.size(), rhs.m_impl.size()));

  auto it = m_impl.begin();
  auto rhs_it = rhs.m_impl.begin();

  while (it != m_impl.end() && rhs_it != rhs.m_impl.end()) {
    if (it->chunk_id < rhs_it->chunk_id) {
      ++it;
      continue;
    }
    if (rhs_it->chunk_id < it->chunk_id) {
      ++rhs_it;
      continue;
    }

    CHECK_EQ(it->chunk_id, rhs_it->chunk_id);
    elem new_elem;
    new_elem.chunk_id = it->chunk_id;
    new_elem.read_id_bits = it->read_id_bits & rhs_it->read_id_bits;
    if (new_elem.read_id_bits) {
      result.m_impl.emplace_back(new_elem);
    }
    ++it;
    ++rhs_it;
  }

  return result;
}

read_id_set& read_id_set::operator|=(const read_id_set& rhs) {
  (*this) = (*this) | rhs;
  return *this;
}

read_id_set& read_id_set::operator&=(const read_id_set& rhs) {
  (*this) = (*this) & rhs;
  return *this;
}

read_id_set& read_id_set::operator-=(const read_id_set& rhs) {
  (*this) = (*this) - rhs;
  return *this;
}

read_id_set read_id_set::operator|(const big_read_id_set& rhs) const {
  read_id_set result = *this;
  result |= rhs;
  return result;
}

read_id_set& read_id_set::operator|=(const big_read_id_set& rhs) {
  big_read_id_set result = rhs;
  result |= *this;
  m_impl.clear();
  m_impl.reserve(result.m_impl.size());
  for (const auto& chunk : result.m_impl) {
    elem new_elem;
    new_elem.chunk_id = chunk.first;
    new_elem.read_id_bits = chunk.second;
    m_impl.push_back(new_elem);
  }
  return *this;
}

read_id_set& read_id_set::operator&=(const big_read_id_set& rhs) {
  auto rhs_it = rhs.m_impl.begin();
  for (auto& chunk : m_impl) {
    auto next_it = rhs.m_impl.find(chunk.chunk_id);
    if (next_it == rhs.m_impl.end()) {
      chunk.read_id_bits = 0;
      continue;
    }

    rhs_it = next_it;
    chunk.read_id_bits &= rhs_it->second;
    ++rhs_it;
  }

  auto new_end = std::remove_if(m_impl.begin(), m_impl.end(),
                                [](const elem& e) -> bool { return e.read_id_bits == 0; });
  m_impl.erase(new_end, m_impl.end());
  return *this;
}

read_id_set read_id_set::operator&(const big_read_id_set& rhs) const {
  read_id_set result;
  result.m_impl.reserve(m_impl.size());

  auto it = m_impl.begin();
  auto rhs_it = rhs.m_impl.begin();

  while (it != m_impl.end() && rhs_it != rhs.m_impl.end()) {
    if (it->chunk_id < rhs_it->first) {
      ++it;
      continue;
    }
    if (rhs_it->first < it->chunk_id) {
      rhs_it = rhs.m_impl.lower_bound(it->chunk_id);
      continue;
    }
    CHECK_EQ(rhs_it->first, it->chunk_id);
    elem new_elem;
    new_elem.read_id_bits = rhs_it->second & it->read_id_bits;
    if (new_elem.read_id_bits) {
      new_elem.chunk_id = rhs_it->first;
      result.m_impl.push_back(new_elem);
    }
    ++it;
    ++rhs_it;
  }
  return result;
}

read_id_set& read_id_set::operator-=(const big_read_id_set& rhs) {
  auto rhs_it = rhs.m_impl.begin();
  for (auto& chunk : m_impl) {
    auto next_it = rhs.m_impl.find(chunk.chunk_id);
    if (next_it == rhs.m_impl.end()) {
      continue;
    }
    rhs_it = next_it;
    chunk.read_id_bits &= ~rhs_it->second;
  }
  auto new_end = std::remove_if(m_impl.begin(), m_impl.end(),
                                [](const elem& e) -> bool { return e.read_id_bits == 0; });
  m_impl.erase(new_end, m_impl.end());
  return *this;
}

read_id_set read_id_set::operator-(const big_read_id_set& rhs) const {
  read_id_set result = *this;
  result -= rhs;
  return result;
}

bool read_id_set::total_order_lt(const read_id_set& rhs) const {
  ssize_t size_diff = ssize_t(m_impl.size()) - ssize_t(rhs.m_impl.size());
  if (size_diff) {
    return size_diff < 0;
  }
  auto it = m_impl.begin();
  auto rhs_it = rhs.m_impl.begin();
  while (it != m_impl.end()) {
    CHECK(rhs_it != rhs.m_impl.end());
    ssize_t chunk_diff = ssize_t(it->chunk_id) - ssize_t(rhs_it->chunk_id);
    if (chunk_diff) {
      return chunk_diff;
    }
    if (it->read_id_bits != rhs_it->read_id_bits) {
      return it->read_id_bits < rhs_it->read_id_bits;
    }

    ++it;
    ++rhs_it;
  }
  CHECK(rhs_it == rhs.m_impl.end());
  return false;
}

read_coverage_t read_coverage_set::build_and_clear(int assembly_len) {
  std::vector<read_coverage_read_t> tot;
  tot.reserve(m_impl.size());
  for (auto& loff_and_reads : m_impl) {
    const auto& loff = loff_and_reads.first;
    read_coverage_read_t new_cov;
    new_cov.read_len = loff.len;
    new_cov.offset = loff.offset;
    new_cov.read_ids = std::move(loff_and_reads.second);
    tot.emplace_back(std::move(new_cov));
  }

  std::sort(tot.begin(), tot.end(), read_coverage_read_order());

  read_coverage_t result(assembly_len, std::move(tot));
  return result;
}

std::vector<int> read_coverage_t::calc_depths(bool include_fwd, bool include_rev, bool interbase,
                                              const readmap* rm) const {
  if (!include_fwd || !include_rev) {
    if (!rm) {
      throw(io_exception("calc_depths: Must specify readmap to do directional calculation"));
    }
  }

  std::vector<int> starts, ends;
  // When calculating, we calculate interbase starts and ends.
  //
  // The first element in starts is how many reads overlap prior to
  // the beginning of the sequence.
  //
  // The second element in starts is how many reads start exactly at
  // the beginning of the sequence.
  //
  // The last element in starts is how many reads start right before
  // the last base of the sequence.
  //
  // The first element in ends is how many reads end right after the
  // first base of the sequence.
  //
  // The second to last element in ends is how many reads end exactly
  // at the end of the sequence.
  //
  // The last element in ends is how many reads overlap past the end
  // of the sequence.
  starts.resize(assembly_len() + 1, 0);
  ends.resize(assembly_len() + 1, 0);

  for (const auto& rd : reads()) {
    CHECK_LT(rd.offset, assembly_len());
    CHECK_GT(rd.offset + rd.read_len, 0);

    auto& starts_here = starts[std::max<int>(rd.offset + 1, 0)];
    auto& ends_here = ends[std::min<int>(rd.offset + rd.read_len - 1, assembly_len())];
    for (uint32_t read_id : rd.read_ids) {
      if (!include_fwd || !include_rev) {
        bool read_is_fwd = rm->get_is_forward(read_id);
        if (read_is_fwd && !include_fwd) {
          continue;
        }
        if (!read_is_fwd && !include_rev) {
          continue;
        }
      }
      ++starts_here;
      ++ends_here;
    }
  }

  std::vector<int> result;

  if (interbase) {
    result.resize(assembly_len() + 1, 0);
    int cur_depth = 0;
    for (unsigned i = 0; i != result.size(); ++i) {
      cur_depth += starts[i];
      result[i] = cur_depth;
      cur_depth -= ends[i];
    }

    CHECK_EQ(0, cur_depth);
  } else {
    result.resize(assembly_len(), 0);

    int cur_depth = starts[0];
    for (unsigned i = 0; i != result.size(); ++i) {
      cur_depth += starts[i + 1];
      result[i] = cur_depth;
      cur_depth -= ends[i];
    }
    cur_depth -= ends[assembly_len()];
    CHECK_EQ(cur_depth, 0);
  }
  return result;
}

read_coverage_t read_coverage_t::get_reads_spanning_offset(int offset) const {
  return get_reads_spanning_offset_internal(offset, false /* don't adjust to 0 */);
}

read_coverage_t read_coverage_t::get_and_adjust_reads_spanning_offset(int offset) const {
  return get_reads_spanning_offset_internal(offset, true /* do adjust to 0 */);
}

void read_coverage_t::adjust_in_place(int offset) {
  for (auto& cov : m_reads) {
    cov.offset += offset;
  }
}

const read_id_set& read_coverage_t::get_read_ids_at(int offset, int read_len) const {
  read_coverage_read_t key;
  key.offset = offset;
  key.read_len = read_len;

  auto it = std::lower_bound(m_reads.begin(), m_reads.end(), key, read_coverage_read_order());
  if (it != m_reads.end() && it->offset == offset && it->read_len == read_len) {
    return it->read_ids;
  }
  static read_id_set empty_read_id_set;
  return empty_read_id_set;
}

read_coverage_t read_coverage_t::subcoverage(int start, size_t len) const {
  int limit = start + len;

  read_coverage_t result;
  result.m_assembly_len = len;

  for (const auto& rd : m_reads) {
    if (rd.offset >= limit) {
      break;
    }
    if (rd.offset + rd.read_len <= start) {
      continue;
    }
    result.m_reads.push_back(rd);
    result.m_reads.back().offset -= start;
  }
  return result;
}

read_coverage_t read_coverage_t::get_reads_spanning_offset_internal(int offset,
                                                                    bool adjust_to_zero) const {
  read_coverage_t result;
  if (!adjust_to_zero) {
    result.m_assembly_len = m_assembly_len;
  }

  result.m_reads.reserve(m_reads.size());
  for (const auto& rd : m_reads) {
    if (rd.offset >= offset) {
      break;
    }
    if (rd.offset + rd.read_len <= offset) {
      continue;
    }
    result.m_reads.push_back(rd);
    if (adjust_to_zero) {
      result.m_reads.back().offset -= offset;
    }
  }
  return result;
}

std::vector<int> read_coverage_t::get_overlaps() const {
  std::vector<int> overlaps;
  overlaps.reserve(m_reads.size());

  if (m_reads.empty()) {
    return overlaps;
  }

  auto it = m_reads.begin();
  int prev_start = it->offset;
  int prev_end = it->offset + it->read_len;

  bool is_first = true;

  for (; it != m_reads.end(); ++it) {
    for (uint32_t read_id __attribute__((unused)) : it->read_ids) {
      if (is_first) {
        is_first = false;
        continue;
      }
      int this_end = it->offset + it->read_len;
      if (prev_start < it->offset && this_end < prev_end) {
        // Subsumed; skip
        continue;
      }

      overlaps.push_back(std::max<int>(prev_end - it->offset, 0));

      prev_start = it->offset;
      prev_end = this_end;
    }
  }

  return overlaps;
}

std::pair<int, int> read_coverage_t::get_overlap_min_max() const {
  int min_ol = std::numeric_limits<int>::max();
  int max_ol = 0;

  if (m_reads.empty()) {
    return std::make_pair(0, 0);
  }

  auto it = m_reads.begin();
  int prev_start = it->offset;
  int prev_end = it->offset + it->read_len;

  bool is_first = true;

  for (; it != m_reads.end(); ++it) {
    // Determine whether there are multiple reads here.
    unsigned read_count = 0;
    for (uint32_t read_id __attribute__((unused)) : it->read_ids) {
      ++read_count;
      if (read_count > 1) {
        // We only care about whether it's 1 or >1 read, so don't count higher than 1.
        break;
      }
    }

    for (size_t i = 0; i < read_count; ++i) {
      if (is_first) {
        is_first = false;
        continue;
      }
      int this_end = it->offset + it->read_len;
      if (prev_start < it->offset && this_end < prev_end) {
        // Subsumed; skip
        continue;
      }

      int overlap = std::max<int>(prev_end - it->offset, 0);

      if (overlap < min_ol) {
        min_ol = overlap;
      }

      if (read_count > 1) {
        // Multiple reads in the same place; they completely overlap each other.
        overlap = it->read_len;

        if (overlap < min_ol) {
          min_ol = overlap;
        }
      }

      if (overlap > max_ol) {
        max_ol = overlap;
      }

      prev_start = it->offset;
      prev_end = this_end;
    }
  }

  if (min_ol == std::numeric_limits<int>::max()) {
    // We didn't find any overlaps.
    min_ol = 0;
  }

  CHECK_GE(max_ol, min_ol);

  return std::make_pair(min_ol, max_ol);
}

size_t read_coverage_t::get_tot_read_count() const {
  size_t tot_count = 0;
  for (const auto& cov_entry : reads()) {
    tot_count += cov_entry.read_ids.size();
  }
  return tot_count;
}

read_id_set read_coverage_t::all_read_ids() const {
  read_id_set tot;

  for (const auto& cov_entry : reads()) {
    tot.insert(cov_entry.read_ids);
  }
  return tot;
}

read_coverage_t read_coverage_t::intersection_with_adjusted(const read_coverage_t& rhs,
                                                            int rhs_adjust) const {
  read_coverage_read_order cov_lt;
  read_coverage_t result;
  result.m_reads.reserve(std::min(m_reads.size(), rhs.m_reads.size()));
  auto it1 = m_reads.begin();
  auto it2 = rhs.m_reads.begin();
  while (it1 != m_reads.end() && it2 != rhs.m_reads.end()) {
    if (cov_lt.lt_with_adjust(*it1, 0, *it2, rhs_adjust)) {
      ++it1;
      continue;
    }

    if (cov_lt.lt_with_adjust(*it2, rhs_adjust, *it1, 0)) {
      ++it2;
      continue;
    }

    CHECK_EQ(it1->offset, it2->offset + rhs_adjust);
    CHECK_EQ(it1->read_len, it2->read_len);

    read_coverage_read_t intersection;
    intersection.read_ids = it1->read_ids.intersection(it2->read_ids);
    if (intersection.read_ids.empty()) {
      ++it1;
      ++it2;
      continue;
    }
    intersection.offset = it1->offset;
    intersection.read_len = it2->read_len;
    result.m_reads.emplace_back(std::move(intersection));
    ++it1;
    ++it2;
  }
  result.m_assembly_len = m_assembly_len;
  return result;
}

read_coverage_t read_coverage_t::union_with(const read_coverage_t& rhs) const {
  read_coverage_read_order cov_lt;
  read_coverage_t result;
  result.m_reads.reserve(m_reads.size() + rhs.m_reads.size());
  auto it1 = m_reads.begin();
  auto it2 = rhs.m_reads.begin();
  while (it1 != m_reads.end() && it2 != rhs.m_reads.end()) {
    if (cov_lt(*it1, *it2)) {
      result.m_reads.push_back(*it1);
      ++it1;
      continue;
    }

    if (cov_lt(*it2, *it1)) {
      result.m_reads.push_back(*it2);
      ++it2;
      continue;
    }

    CHECK_EQ(it1->offset, it2->offset);
    CHECK_EQ(it1->read_len, it2->read_len);

    read_coverage_read_t u;
    u.read_ids = it1->read_ids | it2->read_ids;
    u.offset = it1->offset;
    u.read_len = it2->read_len;
    result.m_reads.emplace_back(std::move(u));
    ++it1;
    ++it2;
  }
  result.m_reads.insert(result.m_reads.end(), it1, m_reads.end());
  result.m_reads.insert(result.m_reads.end(), it2, rhs.m_reads.end());
  if (m_assembly_len == rhs.m_assembly_len) {
    result.m_assembly_len = m_assembly_len;
  }
  return result;
}

read_coverage_t& read_coverage_t::operator|=(const read_coverage_t& rhs) {
  read_coverage_read_order cov_lt;
  read_coverage_t result;
  result.m_reads.reserve(m_reads.size() + rhs.m_reads.size());
  auto it1 = m_reads.begin();
  auto it2 = rhs.m_reads.begin();
  while (it1 != m_reads.end() && it2 != rhs.m_reads.end()) {
    if (cov_lt(*it1, *it2)) {
      result.m_reads.push_back(std::move(*it1));
      ++it1;
      continue;
    }

    if (cov_lt(*it2, *it1)) {
      result.m_reads.push_back(*it2);
      ++it2;
      continue;
    }

    CHECK_EQ(it1->offset, it2->offset);
    CHECK_EQ(it1->read_len, it2->read_len);

    read_coverage_read_t u;
    u.read_ids = it1->read_ids | it2->read_ids;
    u.offset = it1->offset;
    u.read_len = it2->read_len;
    result.m_reads.emplace_back(std::move(u));
    ++it1;
    ++it2;
  }
  result.m_reads.insert(result.m_reads.end(), std::make_move_iterator(it1),
                        std::make_move_iterator(m_reads.end()));
  result.m_reads.insert(result.m_reads.end(), it2, rhs.m_reads.end());
  m_reads = std::move(result.m_reads);
  return *this;
}

read_coverage_t read_coverage_t::difference_with(const read_coverage_t& rhs) const {
  read_coverage_read_order cov_lt;
  read_coverage_t result;
  result.m_reads.reserve(m_reads.size());
  auto it1 = m_reads.begin();
  auto it2 = rhs.m_reads.begin();
  while (it1 != m_reads.end() && it2 != rhs.m_reads.end()) {
    if (cov_lt(*it1, *it2)) {
      result.m_reads.push_back(*it1);
      ++it1;
      continue;
    }

    if (cov_lt(*it2, *it1)) {
      ++it2;
      continue;
    }

    CHECK_EQ(it1->offset, it2->offset);
    CHECK_EQ(it1->read_len, it2->read_len);

    read_coverage_read_t u;
    u.read_ids = it1->read_ids - it2->read_ids;
    u.offset = it1->offset;
    u.read_len = it2->read_len;
    if (!u.read_ids.empty()) {
      result.m_reads.emplace_back(std::move(u));
    }
    ++it1;
    ++it2;
  }
  result.m_reads.insert(result.m_reads.end(), it1, m_reads.end());
  if (m_assembly_len == rhs.m_assembly_len) {
    result.m_assembly_len = m_assembly_len;
  }
  return result;
}

read_coverage_t& read_coverage_t::operator&=(const read_id_set& rhs) {
  for (auto& cov_entry : m_reads) {
    cov_entry.read_ids &= rhs;
  }
  auto new_end = std::remove_if(
      m_reads.begin(), m_reads.end(),
      [](const read_coverage_read_t& cov_entry) { return cov_entry.read_ids.empty(); });
  m_reads.erase(new_end, m_reads.end());
  return *this;
}

read_coverage_t read_coverage_t::operator&(const read_id_set& rhs) const {
  read_coverage_t result = *this;
  result &= rhs;
  return result;
}

read_coverage_t read_coverage_t::operator-(const read_id_set& rhs) const {
  read_coverage_t result = *this;
  result -= rhs;
  return result;
}

read_coverage_t& read_coverage_t::operator-=(const read_id_set& rhs) {
  for (auto& cov_entry : m_reads) {
    cov_entry.read_ids -= rhs;
  }
  auto new_end = std::remove_if(
      m_reads.begin(), m_reads.end(),
      [](const read_coverage_read_t& cov_entry) { return cov_entry.read_ids.empty(); });
  m_reads.erase(new_end, m_reads.end());
  return *this;
}

int read_coverage_t::get_max_flank(int offset) const {
  int max_flank = 0;

  for (const auto& rd : m_reads) {
    if (rd.offset >= offset) {
      break;
    }
    int rd_end = rd.offset + rd.read_len;
    if (rd_end <= offset) {
      continue;
    }

    int left_flank = offset - rd.offset;
    int right_flank = rd_end - offset;
    int flank = std::min(left_flank, right_flank);
    if (flank > max_flank) {
      max_flank = flank;
    }
  }
  return max_flank;
}

std::ostream& operator<<(std::ostream& os, const read_coverage_read_t& rd) {
  return os << "Read@" << rd.offset << "+" << rd.read_len;
}

std::ostream& operator<<(std::ostream& os, const read_coverage_t& rds) {
  const auto& reads = rds.reads();
  os << "[" << reads.size() << " reads: ";
  for (const auto& rd : reads) {
    os << " @" << rd.offset << "+" << rd.read_len;
  }
  os << "]";
  return os;
}

void big_read_id_set::insert(uint32_t read_id) {
  uint32_t chunk_id = read_id / k_mask_bits;
  unsigned offset = read_id & (k_mask_bits - 1);
  m_impl[chunk_id] |= 1ULL << offset;
}

size_t big_read_id_set::size() const {
  size_t tot = 0;
  for (const auto& e : m_impl) {
    tot += __builtin_popcountl(e.second);
  }
  return tot;
}

big_read_id_set& big_read_id_set::operator|=(const read_id_set& rhs) {
  auto it = m_impl.begin();
  for (const auto& chunk : rhs.m_impl) {
    it = m_impl.insert(it, std::make_pair(chunk.chunk_id, 0));
    it->second |= chunk.read_id_bits;
  }
  return *this;
}

big_read_id_set& big_read_id_set::operator-=(const read_id_set& rhs) {
  auto it = m_impl.begin();
  for (const auto& chunk : rhs.m_impl) {
    auto next_it = m_impl.find(chunk.chunk_id);
    if (next_it == m_impl.end()) {
      continue;
    }
    it = next_it;
    it->second &= ~chunk.read_id_bits;
    if (it->second) {
      ++it;
    } else {
      it = m_impl.erase(it);
    }
  }
  return *this;
}

big_read_id_set& big_read_id_set::operator&=(const read_id_set& rhs) {
  auto it = m_impl.begin();
  for (const auto& chunk : rhs.m_impl) {
    auto next_it = m_impl.find(chunk.chunk_id);
    if (next_it == m_impl.end()) {
      continue;
    }
    it = m_impl.erase(it, next_it);
    it->second &= chunk.read_id_bits;
    if (it->second) {
      ++it;
    } else {
      it = m_impl.erase(it);
    }
  }
  m_impl.erase(it, m_impl.end());
  return *this;
}

std::ostream& operator<<(std::ostream& os, const big_read_id_set& ids) {
  bool first = true;
  os << "BigReadIds(";
  for (uint32_t read_id : ids) {
    if (first) {
      first = false;
    } else {
      os << ",";
    }
    os << read_id;
  }
  return os << ")";
}

}  // namespace variants

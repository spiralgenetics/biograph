#include "modules/io/int_map_interface.h"
#include "modules/io/make_unique.h"
#include "modules/io/packed_varbit_vector.h"
#include "modules/io/packed_varint_vector.h"
#include "modules/io/packed_vector.h"
#include "modules/io/parallel.h"

uint8_int_map::uint8_int_map(membuf buffer) : m_buffer(buffer) {}

uint64_t uint8_int_map::get(uint64_t index) const {
  CHECK_LT(index, size());
  return reinterpret_cast<const uint8_t*>(m_buffer.data())[index];
}

uint64_t uint8_int_map::size() const { return m_buffer.size(); }

uint64_t uint8_int_map::max_value() const { return std::numeric_limits<uint8_t>::max(); }

membuf_cachelist uint8_int_map::membufs() const {
  return m_buffer;
}

uint16_int_map::uint16_int_map(membuf buffer) : m_buffer(buffer) {
  CHECK_EQ(0, m_buffer.size() % sizeof(uint16_t));
}

uint64_t uint16_int_map::get(uint64_t index) const {
  CHECK_LT(index, size());
  return reinterpret_cast<const uint16_t*>(m_buffer.data())[index];
}

uint64_t uint16_int_map::size() const { return m_buffer.size() / 2; }

uint64_t uint16_int_map::max_value() const { return std::numeric_limits<uint16_t>::max(); }

membuf_cachelist uint16_int_map::membufs() const {
  return m_buffer;
}

std::unique_ptr<int_map_interface> int_map_interface::detect_subpart(
    const spiral_file_open_state& state) {
  std::vector<std::string> errors;
  try {
    std::unique_ptr<packed_varbit_vector> vb = make_unique<packed_varbit_vector>(state);
    std::unique_ptr<int_map_interface> intf = vb->get_int_map_interface();
    return intf;
  } catch (const io_exception& e) {
    errors.push_back(e.what());
  }
  try {
    return make_unique<packed_varint_vector>(state);
  } catch (const io_exception& e) {
    errors.push_back(e.what());
  }
  try {
    return make_unique<packed_vector<uint32_t, 32>>(state);
  } catch (const io_exception& e) {
    errors.push_back(e.what());
  }

  std::string error_str;
  for (const auto& e : errors) {
    if (!error_str.empty()) {
      error_str += ",";
    }
    error_str += e;
  }
  throw(io_exception("Couldn't audotect int map subpart: " + error_str));
}

std::unique_ptr<int_map_interface> int_map_interface::detect_subpart_or_uint8_membuf(
    const spiral_file_open_state& parent_state, const std::string& subpart_name) {
  if (parent_state.membuf_present(subpart_name)) {
    return make_unique<uint8_int_map>(parent_state.open_membuf(subpart_name));
  }
  return detect_subpart(parent_state.open_subpart(subpart_name));
}

std::unique_ptr<int_map_interface> int_map_interface::detect_subpart_or_uint16_membuf(
    const spiral_file_open_state& parent_state, const std::string& subpart_name) {
  if (parent_state.membuf_present(subpart_name)) {
    return make_unique<uint16_int_map>(parent_state.open_membuf(subpart_name));
  }
  return detect_subpart(parent_state.open_subpart(subpart_name));
}

struct less_than_search::impl_t {
  static constexpr size_t k_factor1 = 64;
  static constexpr size_t k_factor2 = 64;

  impl_t(size_t num_mins, size_t max_value);
  void fill(less_than_search* lt, const int_map_interface* vals);

  mutable_packed_varbit_vector min_vals;
};

constexpr size_t less_than_search::impl_t::k_factor1;
constexpr size_t less_than_search::impl_t::k_factor2;

less_than_search::impl_t::impl_t(size_t num_mins, size_t max_value)
    : min_vals(num_mins, max_value, "less_than_search") {}

void less_than_search::impl_t::fill(less_than_search* lt, const int_map_interface* vals) {
  parallel_for(0, lt->m_num_factor2, [this, lt, vals](size_t start, size_t limit) {
    for (size_t f2_pos = start; f2_pos != limit; ++f2_pos) {
      size_t f1_limit = (f2_pos + 1) * k_factor2;
      if (f1_limit > lt->m_num_factor1) {
        f1_limit = lt->m_num_factor1;
      }
      size_t f2_min = std::numeric_limits<size_t>::max();
      for (size_t f1_pos = f2_pos * k_factor2; f1_pos != f1_limit; ++f1_pos) {
        size_t pos_limit = (f1_pos + 1) * k_factor1;
        if (pos_limit > lt->m_size) {
          pos_limit = lt->m_size;
        }
        size_t f1_min = std::numeric_limits<size_t>::max();
        for (size_t pos = f1_pos * k_factor1; pos != pos_limit; ++pos) {
          f1_min = std::min(f1_min, vals->get(pos));
        }
        min_vals.set(f1_pos, f1_min);
        f2_min = std::min(f1_min, f2_min);
      }
      min_vals.set(lt->m_num_factor1 + f2_pos, f2_min);
    }
  });
}

size_t less_than_search::next_forward_lt(size_t start_pos, size_t max_val) const {
  CHECK_LT(start_pos, m_size);
  size_t pos = start_pos;
  if (m_vals->get(pos) < max_val) {
    return pos;
  }
  while (pos < m_size) {
    ++pos;
    if (pos == m_size) {
      return pos;
    }
    if ((pos % impl_t::k_factor1) == 0) {
      size_t f1_pos = pos / impl_t::k_factor1;
      if ((f1_pos % impl_t::k_factor2) == 0) {
        size_t f2_pos = f1_pos / impl_t::k_factor2;

        if (get_f2(f2_pos) >= max_val) {
          pos = (f2_pos + 1) * impl_t::k_factor1 * impl_t::k_factor2 - 1;
          if (pos >= m_size) {
            return m_size;
          }
          continue;
        }
      }

      if (get_f1(f1_pos) >= max_val) {
        pos = (f1_pos + 1) * impl_t::k_factor1 - 1;
        if (pos >= m_size) {
          return m_size;
        }
        continue;
      }
    }
    if (m_vals->get(pos) < max_val) {
      return pos;
    }
  }
  return pos;
}

size_t less_than_search::next_backward_lt(size_t start_pos, size_t max_val) const {
  CHECK_LT(start_pos, m_size);
  size_t pos = start_pos;
  if (m_vals->get(pos) < max_val) {
    return pos;
  }
  while (pos > 0) {
    --pos;
    if ((pos % impl_t::k_factor1) == (impl_t::k_factor1 - 1)) {
      size_t f1_pos = pos / impl_t::k_factor1;
      if ((f1_pos % impl_t::k_factor2) == (impl_t::k_factor2 - 1)) {
        size_t f2_pos = f1_pos / impl_t::k_factor2;

        if (get_f2(f2_pos) >= max_val) {
          pos = f2_pos * impl_t::k_factor1 * impl_t::k_factor2;
          continue;
        }
      }

      if (get_f1(f1_pos) >= max_val) {
        pos = f1_pos * impl_t::k_factor1;
        continue;
      }
    }
    if (m_vals->get(pos) < max_val) {
      return pos;
    }
  }
  return pos;
}

less_than_search::less_than_search(const int_map_interface* vals)
    : m_vals(vals),
      m_size(vals->size()),
      m_num_factor1((m_size + (impl_t::k_factor1 - 1)) / impl_t::k_factor1),
      m_num_factor2((m_num_factor1 + (impl_t::k_factor2 - 1)) / impl_t::k_factor2),
      m_impl(make_unique<impl_t>(m_num_factor1 + m_num_factor2, vals->max_value())) {
  m_impl->fill(this, vals);
}

less_than_search::~less_than_search() {}

size_t less_than_search::get_f1(size_t pos) const {
  CHECK_LT(pos, m_num_factor1);
  return m_impl->min_vals.get(pos);
}

size_t less_than_search::get_f2(size_t pos) const {
  CHECK_LT(pos, m_num_factor2);
  return m_impl->min_vals.get(pos + m_num_factor1);
}

size_t less_than_search::get_factor1() { return impl_t::k_factor1; }
size_t less_than_search::get_factor2() { return impl_t::k_factor2; }

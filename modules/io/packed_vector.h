#pragma once

// packed_vector provides a mmappable vector containing fixed bit
// values.  Values can contain a number of bits that's any power of 2
// between 1 and 64.

#include "base/base.h"
#include "modules/io/log.h"
#include "modules/io/membuf.h"
#include "modules/io/spiral_file.h"
#include "modules/io/version.h"
#include "modules/io/int_map_interface.h"

#include <string.h>
#include <thread>
#include <vector>

namespace detail {
template <class T>
class accessor {
 public:
  typedef typename T::element_type element_type;
  typedef typename T::value_type value_type;

  accessor(const element_type* buf, size_t pos)
      : m_buf(buf),
        m_bit_offset(pos * T::value_width),
        m_index(m_bit_offset / T::element_width),
        m_start(m_bit_offset % T::element_width),
        m_mask((1ULL << T::value_width) - 1) {}

  operator typename T::value_type() const ATTRIBUTE_NO_SANITIZE_THREAD {
    auto element = m_buf[m_index];
    auto raw = (element >> m_start) & m_mask;
    return typename T::value_type(raw);
  }

  void prefetch_write() const {
    __builtin_prefetch(m_buf + m_index, 1 /* write */);
  }
  void prefetch_read() const {
    __builtin_prefetch(m_buf + m_index, 0 /* read */);
  }

 protected:
  const element_type* const m_buf;
  const size_t m_bit_offset;
  const size_t m_index;
  const size_t m_start;
  const element_type m_mask;
};

template <class T>
class mutable_accessor : public accessor<T> {
 public:
  typedef typename T::element_type element_type;
  typedef typename T::value_type value_type;

  mutable_accessor(element_type* buf, size_t pos)
      : accessor<T>(buf, pos), m_mutable_buf(buf) {}

  void operator=(const typename T::value_type& value) const
      ATTRIBUTE_NO_SANITIZE_THREAD {
    typedef typename T::element_type element_type;
    auto native_value = static_cast<element_type>(value);
    auto append = (native_value << this->m_start);
    auto probe = &m_mutable_buf[this->m_index];
    while (true) {
      volatile auto old_value = *probe;
      auto new_value = (old_value & ~(this->m_mask << this->m_start)) | append;
      if (__sync_bool_compare_and_swap(probe, old_value, new_value)) {
        return;
      }
      std::this_thread::yield();
    }
  }

  // Performs better, but caller guarantees no data contention.
  void set_unlocked(const typename T::value_type& value) const {
    typedef typename T::element_type element_type;
    auto native_value = static_cast<element_type>(value);
    auto append = (native_value << this->m_start);
    auto& elem = m_mutable_buf[this->m_index];
    elem = (elem & ~(this->m_mask << this->m_start)) | append;
  }

  bool compare_and_swap(const typename T::value_type& old_value,
                        const typename T::value_type& new_value) const
      ATTRIBUTE_NO_SANITIZE_THREAD {
    typedef typename T::value_type value_type;
    typedef typename T::element_type element_type;
    element_type* element_ptr = &m_mutable_buf[this->m_index];
    element_type old_element = *element_ptr;
    value_type stored_old_value =
        ((old_element >> this->m_start) & this->m_mask);

    if (stored_old_value != old_value) {
      return false;
    }

    element_type new_element =
        (old_element & ~(this->m_mask << this->m_start)) |
        (element_type(new_value) << this->m_start);
    return __sync_bool_compare_and_swap(element_ptr, old_element, new_element);
  }

  bool safe_increment() const ATTRIBUTE_NO_SANITIZE_THREAD {
    auto probe = &m_mutable_buf[this->m_index];
    while (true) {
      volatile auto old_value = *probe;
      auto raw = (old_value >> this->m_start) & this->m_mask;
      if (raw == T::max_value_static()) {
        return true;
      }
      auto increment = 1ULL << this->m_start;
      auto new_value = old_value + increment;
      if (__sync_bool_compare_and_swap(probe, old_value, new_value)) {
        return false;
      }
    }
  }

 private:
  element_type* const m_mutable_buf;
};

}  // namespace detail

template <typename T, size_t W>
class packed_vector : public int_map_interface {
 public:
  typedef T value_type;
  static const constexpr size_t value_width =
      W;  // size of packed element in bits
  typedef uint64_t element_type;
  static const constexpr size_t element_width = sizeof(element_type) * 8;
  static constexpr size_t values_per_element = element_width / value_width;
  typedef const packed_vector<T, W> this_type;

  static_assert((value_width & (value_width - 1)) == 0,
                "Packed vector currently only supports bit widths that are a "
                "power of 2.");

  typedef const detail::accessor<this_type> const_accessor;

  size_t capacity() const {
    return (m_membuf.size() / sizeof(element_type)) * values_per_element;
  }

  const_accessor at(size_t pos) const {
    return const_accessor(container(), pos);
  }

  const_accessor operator[](size_t pos) const {
    return const_accessor(container(), pos);
  }

  size_t size() const override { return m_size; }
  uint64_t get(uint64_t index) const override {
    return at(index);
  }
  uint64_t max_value() const override {
    return packed_vector::max_value_static();
  }

  size_t memory_used() const { return m_membuf.size(); }

  static constexpr size_t max_value_static() {
    return (value_width == 64) ? std::numeric_limits<uint64_t>::max() : ((1ULL << value_width) - 1);
  }

  static constexpr size_t memory_usage(size_t capacity) {
    return required_elements(capacity) * sizeof(element_type);
  }

  static const product_version packed_vector_version;

  // Opens from a spiral file.  If copy_to_ram is true, creates a copy
  // of this packed vector in RAM.  This can allow for better
  // performance by allowing it to use huge pages.
  packed_vector(const spiral_file_open_state& state,
                std::string description = "(unused)") {
    state.enforce_max_version("packed_vector", packed_vector_version);
    m_membuf = state.open_membuf("packed_data");
    pv_metadata md = state.open_json<pv_metadata>("packed_vector.json");
    m_size = md.value_count;
    CHECK_GT(m_size, 0);
    CHECK_EQ(value_width, md.value_width_bits);
  }
  virtual ~packed_vector() = default;

  membuf_cachelist membufs() const override { return m_membuf; }

 protected:
  // Metadata for this packed vector for when serialized.
  struct pv_metadata {
    TRANSFER_OBJECT {
      VERSION(0);
      FIELD(value_count, TF_STRICT);
      FIELD(value_width_bits, TF_STRICT);
    }

    // Number of values present.
    size_t value_count = 0;

    // Number of bits per value.
    size_t value_width_bits = 0;
  };

  packed_vector() {}

  static constexpr size_t required_elements(size_t capacity) {
    return (value_width * capacity + (element_width - 1)) / element_width;
  }

  const element_type* container() const {
    return reinterpret_cast<const element_type*>(m_membuf.data());
  }

  membuf m_membuf;
  size_t m_size = 0;
};

template <typename T, size_t W>
const product_version packed_vector<T, W>::packed_vector_version{"1.0.0"};

template <typename T, size_t W>
class mutable_packed_vector : public packed_vector<T, W> {
 public:
  typedef mutable_packed_vector<T, W> this_type;
  typedef typename packed_vector<T, W>::element_type element_type;
  typedef const detail::mutable_accessor<this_type> mutable_accessor;

  explicit mutable_packed_vector(size_t size, const std::string& description) {
    this->m_mutable_membuf = new owned_membuf(this->memory_usage(size),
                                              "packed_vector " + description);
    this->m_membuf = m_mutable_membuf;
    this->m_size = size;
    this->m_mutable_membuf.populate_pages_for_write();
  }

  mutable_packed_vector(const spiral_file_open_state& state)
      : packed_vector<T, W>(state) {
    m_mutable_membuf = state.open_mutable_membuf("packed_data");
  }

  mutable_packed_vector(const spiral_file_create_state& state, size_t size) {
    state.set_version("packed_vector", this->packed_vector_version);
    this->m_membuf = m_mutable_membuf =
        state.create_membuf("packed_data", this->memory_usage(size));
    this->m_size = size;
    typename this_type::pv_metadata md;
    md.value_count = this->m_size;
    md.value_width_bits = this->value_width;
    state.create_json<typename this_type::pv_metadata>("packed_vector.json",
                                                       md);
  }

  void reset() {
    memset(m_mutable_membuf.mutable_data(), 0, m_mutable_membuf.size());
  }

  // Equivalent to:
  // while (pos != size() && at(pos).safe_increment()) ++pos;
  size_t claim_next_available(size_t pos) {
    while (pos < this->size()) {
      if (!at(pos).safe_increment()) {
        return pos;
      }
      ++pos;
      size_t bit_offset = pos * this_type::value_width;
      size_t index = bit_offset / this_type::element_width;
      while (this->mutable_container()[index] == ~element_type(0)) {
        // No available free values in this element; skip it entirely.
        ++index;
        bit_offset = index * this_type::element_width;
        pos = bit_offset / this_type::value_width;
        if (pos >= this->size()) {
          return this->size();
        }
      }
    }
    return pos;
  }

  using packed_vector<T, W>::operator[];
  mutable_accessor operator[](size_t pos) {
    return mutable_accessor(mutable_container(), pos);
  }

  using packed_vector<T, W>::at;
  mutable_accessor at(size_t pos) {
    return mutable_accessor(mutable_container(), pos);
  }

 private:
  element_type* mutable_container() const {
    return reinterpret_cast<element_type*>(m_mutable_membuf.mutable_data());
  }
  mutable_membuf m_mutable_membuf;
};

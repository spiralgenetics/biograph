#pragma once

#include "modules/io/int_map_interface.h"
#include "modules/io/membuf.h"
#include "modules/io/spiral_file.h"
#include "modules/io/transfer_object.h"
#include "modules/io/version.h"

#include <mutex>

struct packed_varint_vector_metadata {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(element_count);
    FIELD(max_value);
    FIELD(bytes_per_element);
  }

  size_t element_count = 0;
  uint64_t max_value = 0;
  unsigned bytes_per_element = 0;
};

// packed_varint_vector provides a mmappable vector containing variable length
// integers.
//
// The max value must be specified upon creation, and no values may be
// added that are >= the max value.
class packed_varint_vector : public int_map_interface {
 public:
  packed_varint_vector(const spiral_file_open_state& state);
  virtual ~packed_varint_vector() = default;

  uint64_t size() const override { return m_metadata.element_count; }
  uint64_t max_value() const override { return m_metadata.max_value; }

  uint64_t get(uint64_t index) const override {
    CHECK_LT(index, m_metadata.element_count);

    size_t base_offset = index * m_metadata.bytes_per_element;
    uint64_t value = 0;
    for (unsigned int i = 0; i < m_metadata.bytes_per_element; i++) {
      value |= uint64_t(((const uint8_t*)m_elements.data())[base_offset + i]) << (i * 8);
    }
    return value;
  }

  static unsigned bytes_for_max_value(uint64_t max_value);

  membuf_cachelist membufs() const override { return m_elements; }

 protected:
  packed_varint_vector() = default;

  static const product_version k_varint_vector_version;

  packed_varint_vector_metadata m_metadata;
  membuf m_elements;
};

class mutable_packed_varint_vector : public packed_varint_vector {
 public:
  mutable_packed_varint_vector(const spiral_file_create_state& state, size_t element_count,
                               uint64_t max_value);
  mutable_packed_varint_vector(size_t element_count, uint64_t max_value);

  void set(size_t index, uint64_t value) {
    CHECK_LT(index, m_metadata.element_count);
    CHECK_LE(value, m_metadata.max_value);

    size_t base_offset = index * m_metadata.bytes_per_element;
    for (unsigned int i = 0; i < m_metadata.bytes_per_element; i++) {
      ((uint8_t*)m_mutable_elements.mutable_data())[base_offset + i] = (value >> (i * 8)) & 0xFF;
    }
  }

  bool compare_and_swap(size_t index, uint64_t old_value, uint64_t new_value) {
    CHECK_LT(index, m_metadata.element_count);
    CHECK_LE(new_value, m_metadata.max_value);

    if (get(index) != old_value) return false;
    // TODO(nils): Need more efficient algorithm than a global lock.
    static std::mutex mu;
    std::lock_guard<std::mutex> l(mu);
    if (get(index) != old_value) return false;
    set(index, new_value);
    return true;
  }

 private:
  mutable_membuf m_mutable_elements;
};

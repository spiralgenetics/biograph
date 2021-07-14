#pragma once

#include "modules/io/int_map_interface.h"
#include "modules/io/membuf.h"
#include "modules/io/spiral_file.h"
#include "modules/io/transfer_object.h"
#include "modules/io/version.h"

#include <mutex>

struct packed_varbit_vector_metadata {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(element_count);
    FIELD(max_value);
    FIELD(bits_per_value);
  }

  size_t element_count = 0;
  size_t max_value = 0;
  unsigned bits_per_value = 0;
};

// packed_varbit_vector provides a mmappable vector containing
// variable length integers.  This is similar to packed_varint_vector,
// but values are bit-packed instead of byte-packed.
//
// The max value must be specified upon creation, and no values may be
// added that are >= the max value.
class packed_varbit_vector {
 public:
  static constexpr unsigned k_bits_per_element = sizeof(uint64_t) * 8;

  packed_varbit_vector(const spiral_file_open_state& state);
  virtual ~packed_varbit_vector() = default;

  size_t size() const { return m_metadata.element_count; }
  uint64_t max_value() const { return m_metadata.max_value; }
  uint64_t get(uint64_t index) const {
    CHECK_LT(index, size());
    return m_impl->get(index);
  }

  membuf_cachelist membufs() const { return m_elements; }

  static unsigned bits_for_value(uint64_t max_value);
  static size_t elements_for_values(size_t element_count, unsigned bits_per_value);
  static size_t calc_size(size_t element_count, uint64_t max_value) {
    return sizeof(uint64_t) * elements_for_values(element_count, bits_for_value(max_value));
  }

  class impl_base : public int_map_interface {
   public:
    ~impl_base() = default;
    virtual void varbit_set(mutable_membuf& mb, size_t index, uint64_t value) const = 0;
    size_t size() const { return m_size; }
    uint64_t max_value() const { return m_max_value; }

    membuf_cachelist membufs() const override { return m_elements; }

   protected:
    impl_base(membuf elements, size_t size, uint64_t max_value)
        : m_elements(elements), m_size(size), m_max_value(max_value) {}
    membuf m_elements;
    size_t m_size;
    uint64_t m_max_value;
  };
  std::unique_ptr<int_map_interface> get_int_map_interface() const;

 protected:
  packed_varbit_vector() = default;

  void common_init();
  const uint64_t* elements() const { return reinterpret_cast<const uint64_t*>(m_elements.data()); }

  static const product_version k_varbit_vector_version;

  packed_varbit_vector_metadata m_metadata;
  membuf m_elements;
  size_t m_value_mask = 0;
  std::unique_ptr<impl_base> m_impl;
};

// mutable_packed_varbit_vector::set is not atomic, but it is
// guaranteed to be threadsafe as long as no two threads write to the
// same value at the same time.
class mutable_packed_varbit_vector : public packed_varbit_vector {
 public:
  mutable_packed_varbit_vector(const spiral_file_create_state& state, size_t element_count,
                               uint64_t max_value);
  mutable_packed_varbit_vector(size_t element_count, uint64_t max_value,
                               const std::string& description);

  void set(size_t index, uint64_t value) {
    CHECK_LT(index, size());
    m_impl->varbit_set(m_mutable_elements, index, value);
  }

  // Only use this for testing and benchmarking.
  mutable_membuf get_internal_elements() { return m_mutable_elements; }

 private:
  mutable_packed_varbit_vector() = default;

  uint64_t* mutable_elements() const {
    return reinterpret_cast<uint64_t*>(m_mutable_elements.mutable_data());
  }

  mutable_membuf m_mutable_elements;
};

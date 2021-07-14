#include "modules/io/packed_varbit_vector.h"
#include "modules/io/make_unique.h"
#include "modules/io/parallel.h"

#include <endian.h>

const product_version packed_varbit_vector::k_varbit_vector_version{"1.0.0"};

constexpr unsigned packed_varbit_vector::k_bits_per_element;

namespace {

class varbit_impl_zero : public packed_varbit_vector::impl_base {
 public:
  varbit_impl_zero(membuf mb, size_t size, uint64_t max_value) : impl_base(mb, size, max_value) {}
  uint64_t get(size_t index) const override { return 0; }
  void varbit_set(mutable_membuf& mb, size_t index, uint64_t value) const override {
    DCHECK_EQ(value, 0);
  }

  membuf_cachelist membufs() const override { return membuf_cachelist(); }
};

template <unsigned bytes_per_value>
class varbit_impl_byte_boundaries : public packed_varbit_vector::impl_base {
 public:
  varbit_impl_byte_boundaries(membuf mb, size_t size, uint64_t max_value)
      : impl_base(mb, size, max_value) {}
  uint64_t get(size_t index) const override {
    size_t byte_index = index * bytes_per_value;
    uint64_t result = 0;
    memcpy(&result, m_elements.data() + byte_index, bytes_per_value);
    return result;
  };
  void varbit_set(mutable_membuf& mb, size_t index, uint64_t value) const override {
    size_t byte_index = index * bytes_per_value;
    memcpy(mb.mutable_data() + byte_index, &value, bytes_per_value);
  };
};

template <unsigned min_bytes_per_value>
class varbit_impl_between_bytes : public packed_varbit_vector::impl_base {
 public:
  using elem_t = typename std::conditional<
      (min_bytes_per_value >= 4), uint64_t,
      typename std::conditional<(min_bytes_per_value >= 2), uint32_t, uint16_t>::type>::type;
  varbit_impl_between_bytes(unsigned bits_per_value, membuf mb, size_t size, uint64_t max_value)
      : impl_base(mb, size, max_value),
        m_bits_per_value(bits_per_value),
        m_value_mask((~uint64_t(0)) >> (sizeof(uint64_t) * 8 - bits_per_value)) {
    CHECK_GT(min_bytes_per_value, 0);
    CHECK_LT(bits_per_value, min_bytes_per_value * 8);
    CHECK_GT(bits_per_value, (min_bytes_per_value - 1) * 8);
  }
  uint64_t get(size_t index) const override {
    size_t start_bit = index * m_bits_per_value;
    size_t start_byte = start_bit / 8;
    unsigned start_bit_in_byte = start_bit % 8;
    if (__builtin_expect(min_bytes_per_value < sizeof(elem_t) &&
                             start_byte + sizeof(elem_t) <= m_elements.size(),
                         1)) {
      // A single elem_t encompasses all of the value so we can do it all at once.
      uint64_t val = *reinterpret_cast<const elem_t*>(m_elements.data() + start_byte);
      return (val >> start_bit_in_byte) & m_value_mask;
    }

    uint64_t value = 0;
    memcpy(&value, m_elements.data() + start_byte, min_bytes_per_value);
    if (start_bit_in_byte) {
      value >>= start_bit_in_byte;
      uint8_t last_byte = m_elements.data()[start_byte + min_bytes_per_value];
      unsigned bits_to_shift = min_bytes_per_value * 8 - start_bit_in_byte;
      if (m_bits_per_value > bits_to_shift) {
        value |= (uint64_t(last_byte) << bits_to_shift);
      }
    }
    return value & m_value_mask;
  };

  void varbit_set(mutable_membuf& mb, size_t index, uint64_t value) const override {
    size_t start_bit = index * m_bits_per_value;
    size_t start_byte = start_bit / 8;
    unsigned start_bit_in_byte = start_bit % 8;
    unsigned bits_left = m_bits_per_value;
    uint8_t* ptr = reinterpret_cast<uint8_t*>(mb.mutable_data() + start_byte);
    for (;;) {
      if (!start_bit_in_byte) {
        // Check and see if we can copy parts in without doing a compare-and-swap.
        if (bits_left >= 32) {
          *reinterpret_cast<uint32_t*>(ptr) = value & 0xFFFFFFFF;
          value >>= 32;
          bits_left -= 32;
          ptr += sizeof(uint32_t);
          continue;
        }
        if (bits_left >= 16) {
          *reinterpret_cast<uint16_t*>(ptr) = value & 0xFFFF;
          value >>= 16;
          bits_left -= 16;
          ptr += sizeof(uint16_t);
          continue;
        }
        if (bits_left >= 8) {
          *ptr = value & 0xFF;
          value >>= 8;
          bits_left -= 8;
          ++ptr;
          continue;
        }
      }

      if (!bits_left) {
        return;
      }

      uint8_t mask = 0xFF;
      if (bits_left < 8) {
        mask >>= 8 - bits_left;
      }
      mask <<= start_bit_in_byte;
      uint8_t val_in_byte = (value << start_bit_in_byte);
      uint8_t byte_mask = ~mask;
      for (;;) {
        uint8_t orig_byte = *ptr;
        uint8_t new_byte = ((*ptr) & byte_mask) | val_in_byte;
        if (__sync_bool_compare_and_swap(ptr, orig_byte, new_byte)) {
          break;
        }
      }
      if (bits_left > (8 - start_bit_in_byte)) {
        bits_left -= (8 - start_bit_in_byte);
        value >>= (8 - start_bit_in_byte);
      } else {
        return;
      }
      ++ptr;
      start_bit_in_byte = 0;
    }
  }

 private:
  const unsigned m_bits_per_value;
  const uint64_t m_value_mask;
};

std::unique_ptr<packed_varbit_vector::impl_base> select_varbit_impl(
    unsigned requested_bits_per_value, membuf mb, size_t size, uint64_t max_value) {
  unsigned min_bytes_per_value = (requested_bits_per_value + 7) / 8;
#define SELECT_IMPL(NBYTES)                                                                     \
  if (min_bytes_per_value == NBYTES) {                                                          \
    if (requested_bits_per_value % 8 == 0) {                                                    \
      CHECK_EQ(NBYTES, requested_bits_per_value / 8);                                           \
      return make_unique<varbit_impl_byte_boundaries<NBYTES>>(mb, size, max_value);             \
    } else {                                                                                    \
      return make_unique<varbit_impl_between_bytes<NBYTES>>(requested_bits_per_value, mb, size, \
                                                            max_value);                         \
    }                                                                                           \
  }
  SELECT_IMPL(1);
  SELECT_IMPL(2);
  SELECT_IMPL(3);
  SELECT_IMPL(4);
  SELECT_IMPL(5);
  SELECT_IMPL(6);
  SELECT_IMPL(7);
  SELECT_IMPL(8);
#undef SELECT_IMPL
  CHECK_EQ(requested_bits_per_value, 0) << "min bytes per value: " << min_bytes_per_value;
  return make_unique<varbit_impl_zero>(mb, size, max_value);
}

}  // namespace

unsigned packed_varbit_vector::bits_for_value(size_t value) {
  unsigned bits = 0;
  while (value) {
    ++bits;
    value >>= 1;
  }
  CHECK_LE(bits, 64);
  return bits;
}

size_t packed_varbit_vector::elements_for_values(size_t element_count, unsigned bits_per_value) {
  size_t tot_bits = element_count * bits_per_value;
  size_t tot_uint64s = (tot_bits + 63) / 64;
  return tot_uint64s;
}

void packed_varbit_vector::common_init() {
  CHECK_EQ(m_metadata.bits_per_value, bits_for_value(m_metadata.max_value));
  CHECK_EQ(m_elements.size(), sizeof(uint64_t) * elements_for_values(m_metadata.element_count,
                                                                     m_metadata.bits_per_value));
  m_value_mask = ((~uint64_t(0)) >> (k_bits_per_element - m_metadata.bits_per_value));

  m_impl = select_varbit_impl(m_metadata.bits_per_value, m_elements, size(), max_value());
}

std::unique_ptr<int_map_interface> packed_varbit_vector::get_int_map_interface() const {
  return select_varbit_impl(m_metadata.bits_per_value, m_elements, size(), max_value());
}

packed_varbit_vector::packed_varbit_vector(const spiral_file_open_state& state) {
  state.enforce_max_version("packed_varbit_vector", k_varbit_vector_version);

  m_metadata = state.open_json<packed_varbit_vector_metadata>("packed_varbit_vector.json");
  m_elements = state.open_membuf("elements");

  common_init();
}

mutable_packed_varbit_vector::mutable_packed_varbit_vector(const spiral_file_create_state& state,
                                                           size_t element_count,
                                                           uint64_t max_value) {
  m_metadata.element_count = element_count;
  m_metadata.max_value = max_value;
  m_metadata.bits_per_value = bits_for_value(m_metadata.max_value);

  state.set_version("packed_varbit_vector", k_varbit_vector_version);

  state.create_json<packed_varbit_vector_metadata>("packed_varbit_vector.json", m_metadata);
  m_elements = m_mutable_elements = state.create_membuf(
      "elements",
      sizeof(uint64_t) * elements_for_values(m_metadata.element_count, m_metadata.bits_per_value),
      state.options().with_delayed_write(true));

  common_init();
}

mutable_packed_varbit_vector::mutable_packed_varbit_vector(size_t element_count, uint64_t max_value,
                                                           const std::string& description) {
  m_metadata.element_count = element_count;
  m_metadata.max_value = max_value;
  m_metadata.bits_per_value = bits_for_value(max_value);

  m_elements = m_mutable_elements = mutable_membuf(new owned_membuf(
      sizeof(uint64_t) * elements_for_values(m_metadata.element_count, m_metadata.bits_per_value),
      description));

  common_init();
}



#include "modules/io/packed_varint_vector.h"

const product_version packed_varint_vector::k_varint_vector_version{"1.0.0"};

packed_varint_vector::packed_varint_vector(
    const spiral_file_open_state& state) {
  state.enforce_max_version("packed_varint_vector", k_varint_vector_version);

  m_metadata = state.open_json<packed_varint_vector_metadata>(
      "packed_varint_vector.json");
  m_elements = state.open_membuf("elements");

  CHECK_EQ(m_metadata.bytes_per_element,
           bytes_for_max_value(m_metadata.max_value));
  CHECK_EQ(m_elements.size(),
           m_metadata.bytes_per_element * (m_metadata.element_count + 1));
}

mutable_packed_varint_vector::mutable_packed_varint_vector(
    const spiral_file_create_state& state, size_t element_count,
    uint64_t max_value) {
  m_metadata.element_count = element_count;
  m_metadata.max_value = max_value;
  m_metadata.bytes_per_element = bytes_for_max_value(max_value);

  state.set_version("packed_varint_vector", k_varint_vector_version);

  state.create_json<packed_varint_vector_metadata>("packed_varint_vector.json", m_metadata);
  // Allocate space for one extra element at the end; this allows us
  // to use atomic operations that don't access memory past the end of the
  // buffer.
  m_elements = m_mutable_elements = state.create_membuf(
      "elements", (m_metadata.element_count + 1) * m_metadata.bytes_per_element);
}

mutable_packed_varint_vector::mutable_packed_varint_vector(size_t element_count,
                                                           uint64_t max_value) {
  m_metadata.element_count = element_count;
  m_metadata.max_value = max_value;
  m_metadata.bytes_per_element = bytes_for_max_value(max_value);

  m_elements = m_mutable_elements = mutable_membuf(new owned_membuf(
      (m_metadata.element_count + 1) * m_metadata.bytes_per_element,
      "packed_varint_vector"));
}

unsigned packed_varint_vector::bytes_for_max_value(uint64_t max_value) {
  unsigned n_bytes = 0;

  do {
    n_bytes++;
    max_value >>= 8;
  } while (max_value != 0);

  return n_bytes;
}

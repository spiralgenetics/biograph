#pragma once

#include "modules/io/membuf.h"
#include "modules/io/spiral_file.h"

// This interface can be implemented by anything which exposes a map
// from uint64_t to uint64_t.  It allows easy backward compatibility
// when opening files.
class int_map_interface {
 public:
  int_map_interface(const int_map_interface&) = delete;
  virtual ~int_map_interface() = default;

  virtual uint64_t get(uint64_t index) const = 0;
  virtual size_t size() const = 0;
  virtual uint64_t max_value() const = 0;

  virtual membuf_cachelist membufs() const = 0;

  static std::unique_ptr<int_map_interface> detect_subpart(const spiral_file_open_state& state);

  static std::unique_ptr<int_map_interface> detect_subpart_or_uint8_membuf(
      const spiral_file_open_state& parent_state, const std::string& subpart_name);
  static std::unique_ptr<int_map_interface> detect_subpart_or_uint16_membuf(
      const spiral_file_open_state& parent_state, const std::string& subpart_name);

 protected:
  int_map_interface() = default;
};

class uint8_int_map : public int_map_interface {
 public:
  uint8_int_map(membuf buffer);
  ~uint8_int_map() = default;
  uint64_t get(uint64_t index) const override;
  uint64_t size() const override;
  uint64_t max_value() const override;

  membuf_cachelist membufs() const override;

 private:
  membuf m_buffer;
};

class uint16_int_map : public int_map_interface {
 public:
  uint16_int_map(membuf buffer);
  ~uint16_int_map() = default;
  uint64_t get(uint64_t index) const override;
  uint64_t size() const override;
  uint64_t max_value() const override;

  membuf_cachelist membufs() const override;

 private:
  membuf m_buffer;
};

class less_than_search {
 public:
  less_than_search(const int_map_interface* vals);
  ~less_than_search();

  size_t next_forward_lt(size_t start_pos, size_t max_val) const;
  size_t next_backward_lt(size_t start_pos, size_t max_val) const;

  // Testing access:
  static size_t get_factor1();
  static size_t get_factor2();

 private:
  struct impl_t;

  size_t get_f1(size_t f1_pos) const;
  size_t get_f2(size_t f2_pos) const;

  const int_map_interface* const m_vals;
  const size_t m_size;
  const size_t m_num_factor1;
  const size_t m_num_factor2;

  std::unique_ptr<impl_t> m_impl;
};

#pragma once

// Sparse_multi provides a sparse mapping table.  For each element A
// in [0, A_max), it provides zero or more indexes B in [0, B_max).
//
// Contrast to "bitcount", which provides similar functionality, but
// provides either zero or one index B in [0, B_max).
//
// Similar to how "bitcount" allows a separate array to associate 1
// piece of data with only some of the elements in the source array,
// "sparse_multi" allows a separate array to associate 1 or more
// pieces of data with only some of the elements in the source array.

#include <boost/iterator/iterator_facade.hpp>

#include "modules/io/bitcount.h"
#include "modules/io/spiral_file.h"
#include "modules/io/version.h"

class sparse_multi {
 public:
  // Loads this sparse multi from a spiral file.
  sparse_multi(const spiral_file_open_state& state);

  // Returns a range [x, y) of dest indexes associated with the
  // provided source index.
  std::pair<uint64_t, uint64_t> lookup(uint64_t source_index) const;
  std::pair<uint64_t, uint64_t> lookup_range(uint64_t source_index_start,
                                             uint64_t source_index_limit) const;
  // When lookup doesn't return an empty range, lookup_lower_bound(x) is
  // the same as lookup(x).first.  Otherwise,
  // lookup_lower_bound(x) <= lookup_lower_bound(x+1).
  uint64_t lookup_lower_bound(uint64_t source_index) const;

  uint64_t reverse_lookup(uint64_t idx) const;

  // Returns true if the given dest index is the first one in a group
  // of 1 or more dests associated with a source.
  bool dest_is_first_in_group(uint64_t dest_index) const { return m_dest_to_mid->get(dest_index); }

  size_t source_elem_count() const { return m_source_to_mid->size(); }
  size_t dest_elem_count() const { return m_dest_to_mid->size(); }

  uint64_t lookup_dest_to_mid(uint64_t dest_index) const;
  uint64_t lookup_mid_to_source(uint64_t mid_index) const;

  static const product_version sparse_multi_version;

  // Makes lookup tables to optimize lookup and reverse_lookup.
  void make_find_count_index() {
    m_source_to_mid->make_find_count_index();
    m_dest_to_mid->make_find_count_index();
  }

  class iterator {
   public:
    iterator(const sparse_multi* sm) : m_sm(sm) {}
    std::pair<uint64_t, std::pair<uint64_t, uint64_t>> operator*() const {
      CHECK(m_sm->m_source_to_mid->get(m_source_index));
      CHECK(m_sm->m_dest_to_mid->get(m_dest_index));

      return std::make_pair(m_source_index, std::make_pair(m_dest_index, m_next_dest_index));
    }

    bool operator==(const iterator& other) const {
      return m_sm == other.m_sm && m_source_index == other.m_source_index &&
             m_dest_index == other.m_dest_index && m_next_dest_index == other.m_next_dest_index;
    }
    bool operator!=(const iterator& other) const { return !(*this == other); }

    iterator& operator++() {
      CHECK_LT(m_source_index, m_sm->source_elem_count());
      do {
        m_source_index++;
        if (m_source_index == m_sm->source_elem_count()) {
          break;
        }
      } while (!m_sm->m_source_to_mid->get(m_source_index));

      m_dest_index = m_next_dest_index;
      calculate_next_dest();
      return *this;
    }

    iterator& seek_to(uint64_t source_index);
    iterator& seek_to_begin() {
      if (m_sm->m_source_to_mid->size() == 0) {
        return seek_to_end();
      }
      seek_to(0);
      return *this;
    }

    iterator& seek_to_end() {
      m_source_index = m_sm->source_elem_count();
      m_dest_index = m_next_dest_index = m_sm->dest_elem_count();
      return *this;
    }

   private:
    void calculate_next_dest() {
      m_next_dest_index = m_dest_index;
      if (m_dest_index < m_sm->dest_elem_count()) {
        do {
          m_next_dest_index++;
        } while (m_next_dest_index < m_sm->dest_elem_count() &&
                 !m_sm->m_dest_to_mid->get(m_next_dest_index));
      }
    }

    const sparse_multi* m_sm;
    uint64_t m_source_index;
    uint64_t m_dest_index;
    uint64_t m_next_dest_index;
  };

  iterator begin() const { return iterator(this).seek_to_begin(); }
  iterator end() const { return iterator(this).seek_to_end(); }

  iterator iterator_at_source(uint64_t source_index) const {
    return iterator(this).seek_to(source_index);
  }

  // Returns a list of membufs to cache if memory caching is requested.
  membuf_cachelist membufs() const;

private:
  friend class sparse_multi_builder;
  sparse_multi(std::unique_ptr<bitcount> source_to_mid, std::unique_ptr<bitcount> dest_to_mid);

  std::unique_ptr<bitcount> m_source_to_mid;
  std::unique_ptr<bitcount> m_dest_to_mid;
};

class sparse_multi_builder {
 public:
  // Starts building a new sparse multi table capable of handling the
  // given number of source and dest items.
  sparse_multi_builder(const spiral_file_create_state& state, uint64_t n_source_elems,
                       uint64_t n_dest_elems);

  // Requests a new destination index for the given source index.
  // When building, source indexes must not decrease.  That is, if you
  // call add(x), you are welcome to subsequently call add(x),
  // add(x+1), add(x+5), or add(x+500), but not add(x-1).
  uint64_t add(uint64_t idx);

  // Builds this sparse_multi from an old-style "readmap" gross/fine
  // lookup table.
  void build_from_old_format(const char* gross_read_buf, const char* fine_read_buf);

  // Finalizes this sparse multi and returns a read-only version.
  std::unique_ptr<sparse_multi> finalize();

  ~sparse_multi_builder() {
    if (m_source_to_mid) finalize();
  }

 private:
  uint64_t m_dest_count = 0;
  uint64_t m_last_source_seen = 0;

  std::unique_ptr<bitcount> m_source_to_mid;
  std::unique_ptr<bitcount> m_dest_to_mid;
};

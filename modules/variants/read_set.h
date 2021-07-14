#pragma once

#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>
#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "base/base.h"
#include "modules/bio_base/readmap.h"

namespace variants {

// read_id_set contains a set of read ids.  Since there are some
// seqset entries that have a high quantity of reads associtaed, when
// we track lists of read ids they can sometimes clump up.  Using a
// mask to keep track of which read ids are present saves a lot of
// memory if a lot of nearby read ids are present in a set.
using read_id_mask_t = uint32_t;
class big_read_id_set;
class read_id_set {
  friend class big_read_id_set;

 public:
  static constexpr size_t k_mask_bits = sizeof(read_id_mask_t) * 8;

 protected:
  struct elem {
    uint32_t chunk_id = std::numeric_limits<uint32_t>::max();
    // Bitmasks of read ids present, starting at chunk_id*64.  When
    // building read coverage, we expect relatively dense regions of
    // read ids all in a clump, so this should be a win for storage
    // space.
    read_id_mask_t read_id_bits = 0;

    bool operator<(uint32_t rhs) const { return chunk_id < rhs; }
    bool operator<(const elem& rhs) const { return chunk_id < rhs.chunk_id; }
    bool operator==(const elem& rhs) const {
      return chunk_id == rhs.chunk_id && read_id_bits == rhs.read_id_bits;
    }
  } __attribute__((packed));
  static constexpr size_t k_num_small_elem = 3;
  using impl_t = boost::container::small_vector<elem, k_num_small_elem>;

 public:
  using value_type = const uint32_t;

  read_id_set() = default;
  read_id_set(const big_read_id_set&);

  template <typename Iterator>
  void insert(Iterator begin, const Iterator& end) {
    while (begin != end) {
      insert(*begin);
      ++begin;
    }
  }
  void insert(uint32_t read_id);
  void insert(const read_id_set& old_read_ids);
  void erase(uint32_t read_id);
  void clear();
  bool empty() const { return m_impl.empty(); }
  bool contains(uint32_t read_id) const;

  read_id_set intersection(const read_id_set& rhs) const;

  // Convert to a regular old vector.  TODO(nils): Figure out how to
  // make gtest recognize this as a real container class without
  // having to resort to this.
  std::vector<uint32_t> to_vector() const;

  // Set union.
  read_id_set operator|(const read_id_set& rhs) const;
  // Set difference.
  read_id_set operator-(const read_id_set& rhs) const;
  // Set intersection.
  read_id_set operator&(const read_id_set& rhs) const;

  read_id_set& operator|=(const read_id_set& rhs);
  read_id_set& operator-=(const read_id_set& rhs);
  read_id_set& operator&=(const read_id_set& rhs);

  read_id_set operator|(const big_read_id_set& rhs) const;
  read_id_set operator-(const big_read_id_set& rhs) const;
  read_id_set operator&(const big_read_id_set& rhs) const;
  read_id_set& operator|=(const big_read_id_set& rhs);
  read_id_set& operator-=(const big_read_id_set& rhs);
  read_id_set& operator&=(const big_read_id_set& rhs);

  // Deprecated
  read_id_set operator+(const read_id_set& rhs) const { return (*this) | rhs; }

  class iterator : public std::iterator<std::forward_iterator_tag, const uint32_t> {
   public:
    bool operator==(const iterator& rhs) const {
      return m_cur == rhs.m_cur && m_read_id == rhs.m_read_id;
    }
    bool operator!=(const iterator& rhs) const { return !(*this == rhs); }
    const uint32_t& operator*() const {
      const_cast<iterator*>(this)->advance_to_current();
      return m_read_id;
    }
    iterator& operator++() {
      advance_to_next();
      return *this;
    }
    iterator operator++(int) {
      iterator old = *this;
      ++*this;
      return old;
    }

    iterator(impl_t::const_iterator new_cur) : m_cur(new_cur) {}

   private:
    void advance_to_current() {
      if (m_read_id == std::numeric_limits<uint32_t>::max()) {
        m_bits_left = m_cur->read_id_bits;
        m_read_id = m_cur->chunk_id * k_mask_bits;
      }
      CHECK(m_bits_left);
      unsigned offset = __builtin_ctzl(m_bits_left);
      m_read_id += offset;
      m_bits_left >>= offset;
      CHECK(m_bits_left & 1);
    }

    void advance_to_next() {
      advance_to_current();
      CHECK(m_bits_left & 1);
      m_bits_left >>= 1;
      ++m_read_id;
      if (!m_bits_left) {
        ++m_cur;
        m_read_id = std::numeric_limits<uint32_t>::max();
      }
    }

    impl_t::const_iterator m_cur;
    uint32_t m_read_id = std::numeric_limits<uint32_t>::max();
    read_id_mask_t m_bits_left = 0;
  };

  iterator begin() const { return iterator(m_impl.begin()); }
  iterator end() const { return iterator(m_impl.end()); }

  bool operator==(const read_id_set& rhs) const { return m_impl == rhs.m_impl; }
  bool operator!=(const read_id_set& rhs) const { return m_impl != rhs.m_impl; }

  using const_iterator = iterator;

  size_t size() const;

  // Return a total ordering on read id sets.
  bool total_order_lt(const read_id_set& rhs) const;
  friend std::ostream& operator<<(std::ostream& os, const read_id_set& ids);

 protected:
  impl_t m_impl;
};

// Variant of read_id_set that allows faster operations at the expense
// of slower copies.
class big_read_id_set {
  friend class read_id_set;
  using impl_t = absl::btree_map<uint32_t /* chunk id */, read_id_mask_t>;
  static constexpr size_t k_mask_bits = read_id_set::k_mask_bits;

 public:
  void insert(uint32_t read_id);
  size_t size() const;
  bool empty() const { return m_impl.empty(); }
  big_read_id_set operator|(const read_id_set& rhs) const {
    big_read_id_set result = *this;
    result |= rhs;
    return result;
  }
  big_read_id_set operator&(const read_id_set& rhs) const {
    big_read_id_set result = *this;
    result &= rhs;
    return result;
  }
  big_read_id_set operator-(const read_id_set& rhs) const {
    big_read_id_set result = *this;
    result -= rhs;
    return result;
  }

  big_read_id_set& operator|=(const read_id_set& rhs);
  big_read_id_set& operator&=(const read_id_set& rhs);
  big_read_id_set& operator-=(const read_id_set& rhs);

  class iterator : public std::iterator<std::forward_iterator_tag, const uint32_t> {
   public:
    bool operator==(const iterator& rhs) const {
      return m_cur == rhs.m_cur && m_read_id == rhs.m_read_id;
    }
    bool operator!=(const iterator& rhs) const { return !(*this == rhs); }
    const uint32_t& operator*() const {
      const_cast<iterator*>(this)->advance_to_current();
      return m_read_id;
    }
    iterator& operator++() {
      advance_to_next();
      return *this;
    }
    iterator operator++(int) {
      iterator old = *this;
      ++*this;
      return old;
    }

    iterator(impl_t::const_iterator new_cur) : m_cur(new_cur) {}

   private:
    void advance_to_current() {
      if (m_read_id == std::numeric_limits<uint32_t>::max()) {
        m_bits_left = m_cur->second;
        m_read_id = m_cur->first * k_mask_bits;
      }
      CHECK(m_bits_left);
      unsigned offset = __builtin_ctzl(m_bits_left);
      m_read_id += offset;
      m_bits_left >>= offset;
      CHECK(m_bits_left & 1);
    }

    void advance_to_next() {
      advance_to_current();
      CHECK(m_bits_left & 1);
      m_bits_left >>= 1;
      ++m_read_id;
      if (!m_bits_left) {
        ++m_cur;
        m_read_id = std::numeric_limits<uint32_t>::max();
      }
    }

    impl_t::const_iterator m_cur;
    uint32_t m_read_id = std::numeric_limits<uint32_t>::max();
    read_id_mask_t m_bits_left = 0;
  };

  iterator begin() const { return iterator(m_impl.begin()); }
  iterator end() const { return iterator(m_impl.end()); }

  using const_iterator = iterator;

  friend std::ostream& operator<<(std::ostream& os, const big_read_id_set& ids);

 private:
  impl_t m_impl;
};

struct read_coverage_read_t {
  int offset = std::numeric_limits<int>::max();
  int read_len = std::numeric_limits<int>::max();
  read_id_set read_ids;

  read_coverage_read_t() = default;
  read_coverage_read_t(int offset_, uint32_t read_id_, int read_len_)
      : offset(offset_), read_len(read_len_) {
    read_ids.insert(read_id_);
  }

  bool operator==(const read_coverage_read_t& rhs) const {
    return read_ids == rhs.read_ids && offset == rhs.offset && read_len == rhs.read_len;
  }
  bool operator!=(const read_coverage_read_t& rhs) const { return !(*this == rhs); }

  friend std::ostream& operator<<(std::ostream& os, const read_coverage_read_t& rd);
};

struct read_coverage_read_order {
  bool operator()(const read_coverage_read_t& lhs, const read_coverage_read_t& rhs) const {
    return lt_with_adjust(lhs, 0, rhs, 0);
  }
  bool lt_with_adjust(const read_coverage_read_t& lhs, int lhs_adjust,
                      const read_coverage_read_t& rhs, int rhs_adjust) const {
    if ((lhs.offset + lhs_adjust) != (rhs.offset + rhs_adjust)) {
      return (lhs.offset + lhs_adjust) < (rhs.offset + rhs_adjust);
    } else if (lhs.read_len != rhs.read_len) {
      return lhs.read_len < rhs.read_len;
    }
    return false;
  }
};

class read_coverage_t {
 public:
  read_coverage_t() = default;
  read_coverage_t(int assembly_len, std::vector<read_coverage_read_t> reads)
      : m_assembly_len(assembly_len), m_reads(std::move(reads)) {}
  explicit read_coverage_t(int assembly_len) : m_assembly_len(assembly_len) {}
  read_coverage_t(const read_coverage_t&) = default;
  read_coverage_t(read_coverage_t&&) = default;
  read_coverage_t& operator=(const read_coverage_t&) = default;
  read_coverage_t& operator=(read_coverage_t&&) = default;

  int assembly_len() const { return m_assembly_len; }
  const std::vector<read_coverage_read_t>& reads() const { return m_reads; }
  size_t size() const { return m_reads.size(); }

  bool operator==(const read_coverage_t& rhs) const {
    return assembly_len() == rhs.assembly_len() && reads() == rhs.reads();
  }
  bool operator!=(const read_coverage_t& rhs) const { return !(*this == rhs); }

  read_coverage_t subcoverage(int start, size_t len) const;

  // Calculates coverage depths of the pileup of these reads.
  //
  // If interbase is true, calculate coverage between bases (e.g. edges).  Otherwise, calculate
  // coverage on bases (like how VCF specifies it).
  //
  // If include_fwd is false, skip reads that are facing forward
  // (where is_original_orientation is true).
  //
  // If include_rev is false, skip reads that are facing backwards
  // (where is_original_orientation - false).
  std::vector<int> calc_depths(bool include_fwd = true, bool include_rev = true,
                               bool interbase = true, const readmap* rm = nullptr) const;

  // Returns all reads aligned at the given specific position.
  const read_id_set& get_read_ids_at(int offset, int read_len) const;

  // Returns all the reads that pass through the given offset.
  read_coverage_t get_reads_spanning_offset(int offset) const;

  // Same as get_reads_spanning_offset, but adjusts the returned
  // coverage so that it is centered at 0 instead of at offset.
  read_coverage_t get_and_adjust_reads_spanning_offset(int offset) const;

  void adjust_in_place(int offset);

  // Returns a list of overlap base counts between all the reads.
  std::vector<int> get_overlaps() const;

  // Returns min(get_overlaps()) and max(get_overlaps()), or (0, 0) if
  // no overlaps available.
  std::pair<int, int> get_overlap_min_max() const;

  // Returns the maximum flank length around the given offset.  For
  // each read that spans across the given offset, its flank length is
  // the minimum of the number of bases in that read to the right of
  // the offset and the number of bases in that read to the left of
  // the offset.
  int get_max_flank(int offset) const;

  size_t get_tot_read_count() const;

  // Returns a set of all read ids present, no matter the alignment.
  read_id_set all_read_ids() const;

  // Set operations on read coverage reads.
  read_coverage_t union_with(const read_coverage_t& rhs) const;
  read_coverage_t intersection_with(const read_coverage_t& rhs) const {
    return intersection_with_adjusted(rhs, 0);
  }
  read_coverage_t intersection_with_adjusted(const read_coverage_t& rhs, int adjusted) const;
  read_coverage_t difference_with(const read_coverage_t& rhs) const;

  friend std::ostream& operator<<(std::ostream& os, const read_coverage_t& rds);

  bool empty() const { return m_reads.empty(); }

  // Set operators versus another read coverage.
  read_coverage_t operator|(const read_coverage_t& rhs) const { return this->union_with(rhs); }
  read_coverage_t operator&(const read_coverage_t& rhs) const {
    return this->intersection_with(rhs);
  }
  read_coverage_t operator-(const read_coverage_t& rhs) const { return this->difference_with(rhs); }

  read_coverage_t& operator|=(const read_coverage_t& rhs);
  read_coverage_t& operator&=(const read_coverage_t& rhs) {
    (*this) = (*this) & rhs;
    return *this;
  }
  read_coverage_t& operator-=(const read_coverage_t& rhs) {
    (*this) = (*this) - rhs;
    return *this;
  }

  // Set operators versus a read id set.
  read_coverage_t operator&(const read_id_set& rhs) const;
  read_coverage_t operator-(const read_id_set& rhs) const;
  read_coverage_t& operator&=(const read_id_set& rhs);
  read_coverage_t& operator-=(const read_id_set& rhs);

 private:
  read_coverage_t get_reads_spanning_offset_internal(int offset, bool adjust_to_zero = false) const;

  int m_assembly_len = 0;
  std::vector<read_coverage_read_t> m_reads;
};

// Used for building read_coverage_t.
class read_coverage_set {
 public:
  read_coverage_set() = default;
  read_coverage_set(const read_coverage_set&) = delete;
  read_coverage_set(read_coverage_set&&) = default;
  read_coverage_set& operator=(const read_coverage_set&) = delete;
  read_coverage_set& operator=(read_coverage_set&&) = default;

  void insert(int offset, uint32_t read_id, int read_len);
  void insert(int offset, const read_id_set& read_ids, int read_len);
  void insert(const read_coverage_read_t& val);

  read_coverage_t build_and_clear(int assembly_len);

 private:
  struct len_and_offset {
    int len;
    int offset;

    bool operator==(const len_and_offset& rhs) const noexcept {
      return len == rhs.len && offset == rhs.offset;
    }
    bool operator<(const len_and_offset& rhs) const noexcept {
      if (offset != rhs.offset) {
        return offset < rhs.offset;
      }
      return len < rhs.len;
    }
  };
  struct len_and_offset_hash {
    size_t operator()(const len_and_offset& val) const noexcept {
      size_t result = boost::hash_value(val.len);
      boost::hash_combine(result, val.offset);
      return result;
    }
  };
  using impl_t = absl::flat_hash_map<len_and_offset, read_id_set, len_and_offset_hash>;

  impl_t m_impl;
};

}  // namespace variants

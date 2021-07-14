#pragma once

#include <memory.h>
#include <stdint.h>
#include <stdlib.h>
#include <boost/iterator/iterator_facade.hpp>
#include "modules/io/progress.h"
#include "modules/io/spiral_file.h"

// A utility class that acts as a bitvector, with special support
// for counting the number of 1 bits from bit 0 to bit n, for any n
// less than size.  Uses .25 bits of overhead for each bit
class bitcount {
 public:
  static const product_version bitcount_version;

  // Compute memory usage based on size
  static size_t compute_size(size_t size);

  size_t size() const { return m_nbits; }
  uint64_t total_bits() const { return m_nbits ? count(m_nbits) : 0; }

  // Associate object with a buffer.   Takes a const buffer, but
  // actually allows it to be modified.  Deprecated.
  bitcount(const void* buf, size_t nbits);

  // Creates a new bitcount to track the given number of bits.
  bitcount(size_t nbits);
  bitcount(const spiral_file_create_state& state, size_t nbits);

  // Reads an existing bitcount.
  bitcount(const spiral_file_open_state& state);

  // Setup bitcount
  void init();  // Must be called before using set if using deprecated bitcount
                // constructor.

  void set(size_t i, bool v);  // Set bit i.
  void set_unlocked(size_t i,
                    bool v);  // Set bit i, but caller must ensure there is no data contention.
  // Returns old value.
  bool atomic_exchange(size_t i, bool v);

  // Finish bit setting, no more changes allowed, return total
  size_t finalize(progress_handler_t = null_progress_handler);

  // Use bitcount
  bool get(size_t i) const;      // get bit i
  // Count the number of true (1) bits < i.
  // It is guaranteed that count(size()) == total_bits().
  size_t count(size_t i) const;

  // Given a count, look up the index that generates it.  Runs in log N time.
  // It is guaranteed that find_count(total_bits()) == size().
  size_t find_count(size_t count) const;

  // Iterators that make things look like a vector of totals, size + 1 big
  class total_iterator
      : public boost::iterator_facade<total_iterator, uint64_t const,
                                      std::random_access_iterator_tag,
                                      uint64_t> {
    friend class bitcount;
    friend class boost::iterator_core_access;

   public:
    total_iterator() : m_bc(NULL), m_offset(0) {}

   private:
    total_iterator(const bitcount* bc, int64_t offset)
        : m_bc(bc), m_offset(offset) {}
    uint64_t dereference() const { return m_bc->count(m_offset); }
    bool equal(const total_iterator& rhs) const {
      return m_bc == rhs.m_bc && rhs.m_offset == m_offset;
    }
    void advance(int64_t diff) { m_offset += diff; }
    void increment() { m_offset++; }
    void decrement() { m_offset--; }
    int64_t distance_to(const total_iterator& rhs) const {
      return rhs.m_offset - m_offset;
    }
    const bitcount* m_bc;
    int64_t m_offset;
  };

  total_iterator begin() const { return total_iterator(this, 0); }
  total_iterator end() const { return total_iterator(this, m_nbits + 1); }

  // Generates an index of counts to indexes, which greatly increases
  // the performance of find_count with a small memory penalty.
  void make_find_count_index();

  // Returns a list of membufs to cache if memory caching is requested.
  membuf_cachelist membufs() const;

 private:
  // Metadata for this bitcount for when serialized.
  struct bc_metadata {
    TRANSFER_OBJECT {
      VERSION(0);
      FIELD(nbits, TF_STRICT);
    };

    // Number of bits present in this bitcount.
    size_t nbits = 0;
  };

  // Storage for the actual bits
  static size_t bits_mem_size(size_t nbits);
  // Storage for subaccum
  static size_t subaccum_mem_size(size_t nbits);
  // Storage for accum, needs room 1 extra bit
  static size_t accum_mem_size(size_t nbits);

  // The number of bits
  size_t m_nbits = 0;

  // True if we're allowed to change bits in this bitcount.
  bool m_mutable = false;

  // The actual bits, in groups of 64, low bit = bit 0
  membuf m_bits;
  mutable_membuf m_mutable_bits;
  const size_t* bits() const {
    return reinterpret_cast<const uint64_t*>(m_bits.data());
  }
  size_t* mutable_bits() {
    return reinterpret_cast<uint64_t*>(m_mutable_bits.mutable_data());
  }

  // Local accumulations
  // Accumulations of each group of 64 bits, as 7 bit numbers 0-64 inclusive, in
  // each byte of 64 bit number
  // IE:
  // 0hhhhhhh0ggggggg0fffffff0eeeeeee0ddddddd0ccccccc0bbbbbbb0aaaaaaa
  // That is, the accumulation of m_bits[0] is in the high bits (hhhhhh) of
  // m_subaccum[0].  The data is always accessed as a 64 bit number even though
  // it can also be view as bytes, so that we don't care about endian.
  // We put the elements high bit first to make it easy to 'shift' the elements
  // we don't care about out of the way
  membuf m_subaccum;
  mutable_membuf m_mutable_subaccum;
  const size_t* subaccum() const {
    return reinterpret_cast<const uint64_t*>(m_subaccum.data());
  }
  size_t* mutable_subaccum() {
    return reinterpret_cast<uint64_t*>(m_mutable_subaccum.mutable_data());
  }

  // Global accumulations
  // For each group of 512 bits, the total of all the bits in that group and
  // all of the groups prior.
  membuf m_accum;
  mutable_membuf m_mutable_accum;
  const size_t* accum() const {
    return reinterpret_cast<const uint64_t*>(m_accum.data());
  }
  size_t* mutable_accum() {
    return reinterpret_cast<uint64_t*>(m_mutable_accum.mutable_data());
  }

  // For the find_count lookup table, have an entry for every
  // 2^k_find_count_count_bits counts.
  //
  // For a bitcount with 50% bits set, a value of 11 narrows the
  // search necessary to 2^4=16 "accum" entries.
  static constexpr unsigned k_find_count_count_bits = 11;

  // When using find_count, the maximum number of "accum" entries to
  // search through with a linear search before giving up and
  // reverting to a binary search.  A value of 32 here means that we
  // will only start reverting to binary search if the density is less
  // than 25%.
  static constexpr unsigned k_find_count_max_linear_search = 32;

  // For the find_count lookup table, have each entry point to the
  // index with resolution 2^k_find_count_index_bits.
  static constexpr unsigned k_find_count_index_bits =
      6 /* 2^6 = 64 bits in each entry in m_bits */
      + 3 /* 2^3 = 8 m_bits entries for each entry of "accum" */;

  // Lookup table for find_count.  Maps from (count <<
  // k_find_count_count_bits) to (index << k_find_count_index_bits).
  // index >> k_find_count_index_bits:
  std::vector<uint32_t> m_count_set_to_index;
};

inline void bitcount::set(size_t i, bool v) {
  atomic_exchange(i, v);
}

inline void bitcount::set_unlocked(size_t i, bool v) {
  DCHECK_LT(i, size());
  uint64_t* elem = &mutable_bits()[i / 64];
  uint64_t bit = uint64_t(1) << (i & 63);
  if (v) {
    *elem |= bit;
  } else {
    *elem &= ~bit;
  }
}

inline bool bitcount::atomic_exchange(size_t i, bool v) {
  DCHECK_LT(i, size());
  uint64_t prev;
  if (v) {
    prev = __sync_fetch_and_or(&mutable_bits()[i / 64], (uint64_t(1) << (i & 63)));
  } else {
    prev = __sync_fetch_and_and(&mutable_bits()[i / 64], ~(uint64_t(1) << (i & 63)));
  }
  return prev & (uint64_t(1) << (i & 63));
}

inline bool bitcount::get(size_t i) const {
  DCHECK_LT(i, size());
  return (bits()[i / 64] & (uint64_t(1) << (i & 63))) != 0;
}

inline size_t bitcount::count(size_t i) const {
#if ADDRESS_SANITIZER
  // bitcount::count accesses uninitialized memory when calling
  // count(size()), but it doesn't actually use that data.  But
  // address sanitizer still complains...
  if (i == size()) {
    if (size() == 0) {
      return 0;
    } else {
      return count(i - 1) + (get(i - 1) ? 1 : 0);
    }
  }
  DCHECK_LT(i, size());
#else
  DCHECK_LE(i, size());
#endif
  // First, count the bits int the final 64 elements
  uint64_t x = bits()[i / 64];
  uint64_t mask = (uint64_t(1) << (i & 63)) - 1;
  x &= mask;
  uint64_t bc =
      __builtin_popcount(x >> 32) + __builtin_popcount(x & 0xffffffff);
  // Next, add up the bytes for all prior 64 elements in this group of 512
  size_t group = i / 512;
  size_t subgroup = (i / 64) & 7;
  // Get the groups subaccum and shift away stuff I don't care about
  uint64_t subaccum = this->subaccum()[group];
  // SPLOG("Subaccum = %016lX, subgroup = %d", subaccum, (int) subgroup);
  // WARNING: Dragons be here
  // Shift away unneeded elements
  // I actually do one less shift than I need and finish by zeroing
  // out the low byte This is because in the worst case (most bits
  // set to 1), the folding does actually overflow, but by making
  // sure the empty space is in the *bottom* of the bits the
  // overflow always happens in the top of the 2 final bytes, and is
  // retained during the final addition.  Bit twiddling is complex.
  subaccum >>= 56 - (8 * subgroup);
  subaccum &= 0xffffffffffffff00l;
  // SPLOG("Subaccum shift = %016lX", subaccum);

  // Fold the 8 possible 7 bit values [0-64] and add, making them 4
  // values of 7 bits [0-128]
  subaccum = (subaccum >> 32) + (subaccum & 0xffffffff);
  // Fold again, down to 2 values of 8/9 bits.  Since the low byte
  // is 0, max of low elment is actually 192, thus it fits in 8
  // bits.  Top might be as large as 256, so 9 bits
  subaccum = (subaccum >> 16) + (subaccum & 0xffff);
  // Do final fold, merging our 9 bit number with our 8 bit one
  subaccum = (subaccum >> 8) + (subaccum & 0xff);
  // Compute final sum
  // SPLOG("accum = %ld, subaccum = %ld, bc = %ld", accum[group], subaccum, bc);
  return accum()[group] + subaccum + bc;
}

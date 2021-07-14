#pragma once

#include <boost/iterator/iterator_facade.hpp>

#include "modules/bio_base/dna_sequence.h"
#include "modules/io/membuf.h"
#include "modules/io/parallel.h"
#include "modules/io/spiral_file.h"
#include "vendor/libdivide/libdivide.h"

namespace build_seqset {

template <typename CounterType>
class kmer_count_table {
 public:
  static constexpr kmer_t k_unused_entry = std::numeric_limits<kmer_t>::max();
  static constexpr kmer_t k_kmer_mask = std::numeric_limits<kmer_t>::max() >> 2;
  static constexpr kmer_t k_fwd_flag = 1ULL << 63;
  static constexpr kmer_t k_rev_flag = 1ULL << 62;

  struct element {
    kmer_t kmer_and_flags = k_unused_entry;
    CounterType fwd_count = 0;
    CounterType rev_count = 0;

    kmer_t kmer() const { return kmer_and_flags & k_kmer_mask; }
    bool fwd_flag() const { return kmer_and_flags & k_fwd_flag; }
    bool rev_flag() const { return kmer_and_flags & k_rev_flag; }
    operator bool() const { return kmer_and_flags != k_unused_entry; }
  };

  static uint64_t hash_kmer(kmer_t kmer) {
    return kmer * 15674341118187572551ULL;
  }

  // "description" is included in the "table-too-small" error to
  // identify which table this is.
  kmer_count_table(size_t table_size, std::string description)
      : m_table_size(std::max(size_t(3), table_size)),
        m_divider(m_table_size),
        m_mutable_table_buffer(new owned_membuf(
            table_size * sizeof(element), "kmer_count_table: " + description)),
        m_mutable_table(
            reinterpret_cast<element*>(m_mutable_table_buffer.mutable_data())),
        m_table_buffer(m_mutable_table_buffer),
        m_table(m_mutable_table),
        m_description(description) {
    for (size_t idx = 0; idx != table_size; ++idx) {
      m_mutable_table[idx].kmer_and_flags = k_unused_entry;
    }
  }

  // Increments the count for the given kmer.  Returns the old value.
  CounterType increment(kmer_t kmer, bool flipped, bool set_fwd_flag = false,
                        bool set_rev_flag = false) {
    CHECK(!m_compacted);
    DCHECK_NE(kmer, k_unused_entry);
    DCHECK_LE(kmer, k_kmer_mask);
    size_t table_pos = modulo_size(hash_kmer(kmer));
    bool wrapped = false;

    auto it = m_mutable_table + table_pos;
    while ((it->kmer_and_flags & k_kmer_mask) != kmer) {
      if (it->kmer_and_flags == k_unused_entry) {
        do {
          __sync_bool_compare_and_swap(&it->kmer_and_flags, k_unused_entry,
                                       kmer);
        } while (it->kmer_and_flags == k_unused_entry);
        continue;
      }

      ++it;
      if (it == (m_mutable_table + m_table_size)) {
        if (wrapped) {
          throw(io_exception("Kmer table (" + m_description + ") too small"));
        }
        wrapped = true;
        it = m_mutable_table;
      }
    }

    if (flipped) {
      std::swap(set_fwd_flag, set_rev_flag);
    }
    kmer_t new_flags =
        (set_fwd_flag ? k_fwd_flag : 0) | (set_rev_flag ? k_rev_flag : 0);

    while (it->kmer_and_flags != (it->kmer_and_flags | new_flags)) {
      kmer_t old_val = it->kmer_and_flags;
      kmer_t new_val = old_val | new_flags;
      __sync_bool_compare_and_swap(&it->kmer_and_flags, old_val, new_val);
    }

    CounterType* counter = (flipped ? &it->rev_count : &it->fwd_count);
    CounterType old_val;
    do {
      old_val = *counter;
      if (old_val == std::numeric_limits<CounterType>::max()) {
        return old_val;
      }
    } while (!__sync_bool_compare_and_swap(counter, old_val, old_val + 1));
    return old_val;
  }

  void prefetch_write(kmer_t kmer) const {
    size_t table_pos = modulo_size(hash_kmer(kmer));
    __builtin_prefetch(m_mutable_table + table_pos, 1 /* write */);
  }

  const element& get(kmer_t kmer) const {
    CHECK(!m_compacted);
    bool wrapped = false;
    size_t table_pos = modulo_size(hash_kmer(kmer));
    auto it = m_table + table_pos;
    while (it->kmer() != kmer && it->kmer_and_flags != k_unused_entry) {
      ++it;
      if (it == (m_table + m_table_size)) {
        if (wrapped) {
          throw(io_exception("Kmer table (" + m_description + ") too small"));
        }
        wrapped = true;
        it = m_table;
      }
    }
    return *it;
  }

  void sort();
  void compact();
  void sort_and_compact() {
    compact();
    sort();
  }

  kmer_count_table(const spiral_file_open_state& state) {
    state.enforce_ephemeral_version("kmer_count_table" +
                                    std::to_string(sizeof(element)));
    m_table_buffer = state.open_membuf("elements");
    m_table = reinterpret_cast<const element*>(m_table_buffer.data());
    CHECK_EQ(0, m_table_buffer.size() % sizeof(element));
    m_table_size = m_table_buffer.size() / sizeof(element);
    m_divider =
        libdivide::divider<uint64_t, libdivide::BRANCHFREE>(m_table_size);

    m_compacted = true;
  }

  void save(const spiral_file_create_state& state) {
    state.set_ephemeral_version("kmer_count_table" +
                                std::to_string(sizeof(element)));
    CHECK(m_compacted);
    state.create_membuf("elements", m_table_buffer);
  }

  using const_iterator = const element*;
  const_iterator begin() const { return m_table; }
  const_iterator end() const { return begin() + m_table_size; }
  size_t size() const { return m_table_size; }

  static constexpr size_t max_value() {
    return std::numeric_limits<CounterType>::max();
  }

 private:
  // Modulo the given hash to [0, size()).  But don't divide; dividing
  // is slow.
  size_t modulo_size(size_t hash) const {
    size_t quotient = m_divider.perform_divide(hash);
    size_t remainder = hash - quotient * m_table_size;
    DCHECK_LT(remainder, m_table_size);
    return remainder;
  }

  bool m_sorted = false;
  bool m_compacted = false;
  size_t m_table_size = 0;

  libdivide::divider<uint64_t, libdivide::BRANCHFREE> m_divider;

  mutable_membuf m_mutable_table_buffer;
  element* m_mutable_table = nullptr;

  membuf m_table_buffer;
  const element* m_table = nullptr;

  std::string m_description;
};

template <typename CounterType>
void kmer_count_table<CounterType>::sort() {
  CHECK(!m_sorted);
  CHECK(m_compacted);
  parallel_sort_in_place(m_mutable_table, m_mutable_table + m_table_size,
                         [](const element& a, const element& b) -> bool {
                           return a.kmer() < b.kmer();
                         });
  m_sorted = true;
}

template <typename CounterType>
void kmer_count_table<CounterType>::compact() {
  CHECK(!m_sorted);
  CHECK(!m_compacted);

  element* orig_end = m_mutable_table + m_table_size;
  element* new_end =
      std::partition(m_mutable_table, orig_end, [](const element& k) -> bool {
        return k.kmer_and_flags != k_unused_entry;
      });
  size_t new_size = new_end - m_mutable_table;
  CHECK_LE(new_size, m_table_size);
  m_table_size = new_size;
  m_mutable_table_buffer.discard_region(
      reinterpret_cast<char*>(new_end),
      reinterpret_cast<char*>(orig_end) - reinterpret_cast<char*>(new_end));
  m_mutable_table_buffer = m_mutable_table_buffer.subbuf(
      0, reinterpret_cast<char*>(new_end) -
             reinterpret_cast<char*>(m_mutable_table));
  m_table = m_mutable_table;
  m_table_buffer = m_mutable_table_buffer;
  CHECK_EQ(m_table_buffer.size(), m_table_size * sizeof(element));
  m_compacted = true;
}

}  // namespace build_seqset

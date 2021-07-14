#pragma once

#include <cstdint>
#include <mutex>
#include "modules/bio_base/dna_sequence.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/string_view.h"
#include "modules/io/parallel.h"

namespace build_seqset {

class seq_repository {
 public:
  static constexpr unsigned k_inline_base_bytes = 7;
  static constexpr unsigned k_offset_and_rc_bytes = 5;

  static constexpr unsigned k_inline_bases = k_inline_base_bytes * 4;
  static constexpr unsigned k_offset_bits = k_offset_and_rc_bytes * 8 - 1;

  static constexpr size_t k_max_offset =
      ~(std::numeric_limits<size_t>::max() << k_offset_bits);

  class entry_data;
  class entry_data_size_checker;

  class entry_base;
  class entry;
  class reference;
  class popped_reference;

  class iterator;
  class popped_iterator;

  class repo_builder;
  class ref_builder;

  seq_repository(const std::string& ref_filename,
                 const std::string& repo_filename);
  seq_repository(const std::string& ref_filename, const dna_slice& repo);
  ~seq_repository();

  const entry_data* data_begin() const;
  const entry_data* data_end() const;
  entry_data* data_begin();
  entry_data* data_end();

  iterator begin() const;
  iterator end() const;
  size_t size() const { return m_size; }

  const dna_slice& repo() const { return m_repo_slice; }

  std::function<bool(const entry_data&, const entry_data&)>
  less_than_using_repo() const;
  std::function<bool(const entry_data&, const dna_slice&)>
  less_than_slice_using_repo() const;

  void sort_entry_data(entry_data* begin, entry_data* end) const;

  void set_delete_on_close(bool new_delete_on_close) {
    m_delete_on_close = new_delete_on_close;
  }

  using value_type = entry;

 private:
  // Provide an ordering on entry data.
  bool compare_entry_data_lt(const entry_data& lhs, const entry_data& rhs) const;

  std::unique_ptr<mmap_buffer> m_ref;
  entry_data* m_entry_data_start = nullptr;
  size_t m_size = 0;

  boost::optional<mmap_buffer> m_repo;
  dna_slice m_repo_slice;

  std::string m_ref_filename;
  bool m_delete_on_close = false;

  static std::atomic<size_t> m_fast_compare_counter;
  static std::atomic<size_t> m_slow_compare_counter;
};

class seq_repository::entry_data {
 public:
  entry_data(unsigned size, dna_slice inline_seq, size_t offset,
             bool offset_is_rc) {
    CHECK_LE(size, std::numeric_limits<decltype(m_size)>::max());
    m_size = size;
    CHECK_EQ(inline_seq.size(), std::min<unsigned>(k_inline_bases, m_size));
    set_inline_bases(inline_seq);
    set_offset_and_rc(offset, offset_is_rc);
  }
  entry_data() : entry_data(0, dna_slice(), k_max_offset, false) {}
  entry_data(const entry_data&) = default;
  entry_data(entry_data&&) = default;
  entry_data& operator=(const entry_data&) = default;
  entry_data& operator=(entry_data&&) = default;

  unsigned size() const { return m_size; }

  dna_slice inline_bases() const {
    return dna_slice(dna_const_iterator(m_inline_bases, 0, 0 /* rev comp */),
                     std::min<unsigned>(k_inline_bases, m_size));
  }

  const uint8_t* raw_inline_bases() const { return m_inline_bases; }
  bool inline_bases_lt(const seq_repository::entry_data& rhs) const {
    return inline_bases_cmp(rhs) < 0;
  }

  int inline_bases_cmp(const seq_repository::entry_data& rhs) const {
    int cmp = memcmp(raw_inline_bases(), rhs.raw_inline_bases(), k_inline_base_bytes);
    if (cmp) {
      // Inline bases differ; sort by them
      return cmp;
    }
    int ptr_cmp = memcmp(m_offset_and_rc, rhs.m_offset_and_rc, k_offset_and_rc_bytes);
    if (ptr_cmp == 0 || size() <= k_inline_bases || rhs.size() <= k_inline_bases) {
      // At least one is inline, or they point to the same place in
      // the sequence repo; we can compare just based on size.
      return int(size()) - int(rhs.size());
    }
    // Otherwise, can't determine order without comparing the bases in the repo.
    return 0;
  }

  std::pair<size_t, bool> offset_and_rc() const {
    size_t offset_and_rc = 0;
    for (uint8_t c : m_offset_and_rc) {
      offset_and_rc <<= 8;
      offset_and_rc |= c;
    }
    return std::make_pair(offset_and_rc >> 1, offset_and_rc & 1);
  }

  // Shifts inline bases to the left (e.g. pop front), adds new_base onto the
  // back.
  void shift_inline_bases(dna_base new_base);

  friend std::ostream& operator<<(std::ostream& os, const entry_data& e) {
    os << "[" << e.size() << " bases, inline=" << e.inline_bases();
    size_t offset;
    bool rc;
    std::tie(offset, rc) = e.offset_and_rc();
    os << " offset=" << offset << (rc ? " REV" : " FWD");
    os << "]";
    return os;
  }

 private:
  void set_inline_bases(dna_slice seq) {
    CHECK_LE(seq.size(), k_inline_bases);
    auto it = seq.begin();
    for (size_t offset = 0; offset < k_inline_base_bytes; offset += sizeof(uint64_t)) {
      uint64_t val = 0;
      unsigned shift = sizeof(uint64_t) * 8;
      for (; it != seq.end() && shift; ++it) {
        shift -= 2;
        val |= uint64_t(int(dna_base(*it))) << shift;
      }
      val = htobe64(val);
      memcpy(m_inline_bases + offset, &val,
             std::min<size_t>(k_inline_base_bytes - offset, sizeof(uint64_t)));
    }
  }

  void set_offset_and_rc(size_t offset, bool rc) {
    size_t offset_and_rc = (offset << 1) + (rc ? 1 : 0);
    for (unsigned i = 0; i < k_offset_and_rc_bytes; ++i) {
      m_offset_and_rc[k_offset_and_rc_bytes - i - 1] = offset_and_rc & 0xFF;
      offset_and_rc >>= 8;
    }
    CHECK_EQ(0, offset_and_rc);
  }

 private:
  friend class seq_repository::entry_data_size_checker;
  uint16_t m_size;
  uint8_t m_inline_bases[k_inline_base_bytes];
  uint8_t m_offset_and_rc[k_offset_and_rc_bytes];
} __attribute__((packed));

class seq_repository::entry_data_size_checker {
  static_assert(sizeof(seq_repository::entry_data) ==
                sizeof(seq_repository::entry_data::m_size) + seq_repository::k_inline_base_bytes +
                seq_repository::k_offset_and_rc_bytes,
                "entry_data should not be padded as this wastes space");
};

class seq_repository::entry_base
    : public dna_sequence_ordered<entry_base, entry_base>,
      public dna_sequence_ordered<entry_base, dna_slice> {
 public:
  dna_sequence sequence() const;
  unsigned size() const { return get_entry_data().size() - popped_count(); }

  dna_compare_result compare_to(const entry_base& rhs) const;

  dna_compare_result compare_to(const dna_slice& rhs) const;

  unsigned shared_prefix_length(const entry_base& rhs) const;

  entry pop_front() const;

  dna_slice get_repo_seq() const;
  virtual const entry_data& get_entry_data() const = 0;
  virtual unsigned popped_count() const { return 0; }
  const dna_slice& get_repo() const { return m_repo; }

  entry_data reify_pop() const;

  friend std::ostream& operator<<(std::ostream& os, const entry_base& eb) {
    os << eb.sequence() << " (@" << eb.get_entry_data().offset_and_rc().first
       << ", pop=" << eb.popped_count() << ")";
    return os;
  }

 protected:
  // Protect all the things.
  entry_base() = default;
  entry_base(const entry_base&) = default;
  entry_base(entry_base&&) = default;
  entry_base& operator=(const entry_base&) = default;
  entry_base& operator=(entry_base&&) = default;
  ~entry_base() = default;

  explicit entry_base(const dna_slice& repo) : m_repo(repo) {}

  void check_same_repo(const entry_base& rhs) const;

  dna_slice m_repo;

 private:
  boost::optional<dna_compare_result> fast_compare_to(const entry_base& rhs) const;
  dna_compare_result slow_compare_to(const entry_base& rhs) const;
};

// Self contained, assignable value type.
class seq_repository::entry final : public entry_base {
 public:
  entry() = default;
  entry(const entry_data& the_entry, const dna_slice& repo, unsigned popped)
      : entry_base(repo), m_data(the_entry), m_popped(popped) {}
  entry(const entry&) = default;
  entry(entry&&) = default;
  entry& operator=(const entry&) = default;
  entry& operator=(entry&&) = default;

  entry(const entry_base&);
  entry& operator=(const entry_base&);

  const entry_data& get_entry_data() const override { return m_data; }
  unsigned popped_count() const override { return m_popped; }

 private:
  entry_data m_data;
  unsigned m_popped = 0;
};

class seq_repository::reference final : public entry_base {
 public:
  reference() = default;
  reference(const entry_data* the_entry, dna_slice repo)
      : entry_base(repo), m_data(the_entry) {}
  reference(const reference&) = default;

  const entry_data& get_entry_data() const override {
    DCHECK(m_data);
    return *m_data;
  }

 protected:
  const entry_data* const m_data = nullptr;
};

class seq_repository::popped_reference final : public entry_base {
 public:
  popped_reference() = default;
  popped_reference(const entry_data* the_entry, dna_slice repo, unsigned popped)
      : entry_base(repo),
        m_data(the_entry),
        m_popped(std::min<unsigned>(popped, m_data->size())) {}
  popped_reference(const popped_reference&) = default;

  const entry_data& get_entry_data() const override {
    DCHECK(m_data);
    return *m_data;
  }
  unsigned popped_count() const override { return m_popped; }

 protected:
  const entry_data* const m_data = nullptr;
  unsigned const m_popped = 0;
};

class seq_repository::popped_iterator
    : public boost::iterator_facade<popped_iterator, entry,
                                    std::random_access_iterator_tag,
                                    popped_reference, ptrdiff_t> {
 public:
  popped_iterator() = default;
  popped_iterator(const popped_iterator&) = default;
  popped_iterator& operator=(const popped_iterator&) = default;
  popped_iterator(const entry_data* pos, const dna_slice& repo, unsigned popped)
      : m_pos(pos), m_repo(repo), m_popped(popped) {}
  popped_iterator pop_front() {
    return popped_iterator(m_pos, m_repo, m_popped + 1);
  }

 private:
  friend class boost::iterator_core_access;

  popped_reference dereference() const {
    DCHECK(m_pos);
    return popped_reference(m_pos, m_repo, m_popped);
  }
  bool equal(const popped_iterator& rhs) const {
    CHECK_EQ(m_popped, rhs.m_popped);
    return rhs.m_pos == m_pos;
  }
  void increment() {
    DCHECK(m_pos);
    m_pos++;
  }
  void decrement() {
    DCHECK(m_pos);
    m_pos--;
  }
  void advance(ptrdiff_t n) {
    DCHECK(m_pos);
    m_pos += n;
  }
  ptrdiff_t distance_to(const popped_iterator& rhs) const {
    DCHECK(m_pos);
    DCHECK(rhs.m_pos);

    return rhs.m_pos - m_pos;
  }

  const entry_data* m_pos = nullptr;
  dna_slice m_repo;
  unsigned m_popped = 0;
};

class seq_repository::iterator
    : public boost::iterator_facade<iterator, entry,
                                    std::random_access_iterator_tag, reference,
                                    ptrdiff_t> {
 public:
  iterator() = default;
  iterator(const iterator&) = default;
  iterator& operator=(const iterator&) = default;

  iterator(const entry_data* pos, dna_slice repo) : m_pos(pos), m_repo(repo) {}

  popped_iterator pop_front() const {
    return popped_iterator(m_pos, m_repo, 1);
  }

 private:
  friend class boost::iterator_core_access;

  reference dereference() const {
    DCHECK(m_pos);
    return reference(m_pos, m_repo);
  }
  bool equal(const iterator& rhs) const { return rhs.m_pos == m_pos; }
  void increment() {
    DCHECK(m_pos);
    m_pos++;
  }
  void decrement() {
    DCHECK(m_pos);
    m_pos--;
  }
  void advance(ptrdiff_t n) {
    DCHECK(m_pos);
    m_pos += n;
  }
  ptrdiff_t distance_to(const iterator& rhs) const { return rhs.m_pos - m_pos; }

  const entry_data* m_pos = nullptr;
  dna_slice m_repo;
};

class seq_repository::repo_builder {
 public:
  repo_builder(const std::string& filename);
  ~repo_builder();

  // Returns the offset of the written sequence, in bases.
  size_t write_seq(const dna_slice& seq);

 private:
  void write_seq_unlocked(dna_slice seq);
  void flush_write_buffer(bool final_flush = false);

  static constexpr size_t k_write_buffer_size = 65536;

  std::mutex m_mu;
  file_writer m_writer;

  size_t m_cur_offset = 0;

  std::unique_ptr<unsigned char[]> m_write_buffer;
  dna_sequence::iterator m_write_buffer_it;
  size_t m_bases_avail = 0;
};

class part_counts;
class seq_repository::ref_builder {
 public:
  ref_builder(const std::string& filename, part_counts* counts = nullptr);
  ~ref_builder();

  void write_sequence(dna_slice seq, seq_repository::repo_builder& repo,
                      unsigned fwd_suffixes = 1, unsigned rc_suffixes = 0);
  void write_entry(const entry_data& the_entry);
  // If do_force is false, write_entries_and_clear will only write the
  // entries if it can easily get the lock.  If do_force is true, it
  // will block until the lock is available.
  //
  // If the entries are written, entries is cleared.
  void write_entries_and_clear(std::vector<entry_data>& entries, bool do_force);
  void write_entry(const entry_base& the_entry);
  // With write_*_unlocked, caller must ensure that multiple
  // threads don't write to the same ref builder at the same time.
  void write_entry_unlocked(const entry_base& the_entry);
  void write_entry_unlocked(const entry_data& the_entry);
  void flush();

 private:
  void write_entries_unlocked(const std::vector<entry_data>& entries);

  std::mutex m_mu;
  file_writer m_writer;
  part_counts* m_part_counts = nullptr;

  static constexpr size_t k_write_buffer_entries = 4096;
  std::vector<entry_data> m_write_buffer;
};

}  // namespace build_seqset

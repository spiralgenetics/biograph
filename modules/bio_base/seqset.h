
#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/seqset_bitmap.h"
#include "modules/io/bitcount.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/make_unique.h"
#include "modules/io/spiral_file.h"
#include "modules/io/version.h"
#include "modules/io/progress.h"
#include "modules/io/packed_varbit_vector.h"
#include "modules/io/int_map_interface.h"
#include <map>
#include <boost/optional.hpp>

class seqset_range;
class readmap;

class seqset {
  friend class seqset_range;

  static const product_version seqset_version;

public:
  static size_t compute_old_format_size(size_t entries);

  ~seqset();               // Clears any internal storage, doesn't free buffer

  // Create a new spiral file
  seqset(const spiral_file_create_state &state, size_t entries, unsigned max_entry_len);
  // Load from an existing spiral file
  seqset(const spiral_file_open_state &state);
  // Autodetect old non-spiral-file format or spiral file:
  seqset(const std::string& path, const spiral_file_options& options);
  seqset(const std::string& path) : seqset(path, spiral_file_options()) {}

  const std::string& path() const { return m_path; }

  // Compatibility with old deprecated seqset_file interface.  Deprecated.
  const seqset& get_seqset() const { return *this; }

  seqset(const seqset &) = delete;
  seqset &operator=(const seqset &) = delete;
  seqset(seqset &&) = delete;
  seqset &operator=(seqset &&) = delete;

  // The following are only used for seqset creation.  TODO(nils):
  // Move these to a separate class.
  void init(); // Clear a pwbt for initial construction
  void set_shared(size_t row, unsigned shared);
  void set_entry_size(size_t row, unsigned entry_size);
  void set_bit(size_t row, dna_base b, bool is_set);
  bitcount& mutable_prev(dna_base b);

  void finalize(
      progress_handler_t prog = null_progress_handler); // Complete the seqset
  // End of seqset creation calls.

  // Returns a list of membufs to cache if memory caching is requested.
  membuf_cachelist membufs() const;

  size_t size() const { return m_entries; }

  // TODO(nils): Remove references to read_len since they won't work
  // with variable length reads.
  // (Adam) It'd still be useful to have read_len mean 'maximum read length'
  // (nils) Also, 'minimum read length'.
  unsigned read_len() const { return m_read_len; }

  unsigned max_read_len() const { return m_entry_sizes->max_value(); }

  // Returns the first range of size 'size', returns invalid if none exists
  seqset_range ctx_begin() const;
  // Makes a range for a specific entry of the table
  inline seqset_range ctx_entry(uint64_t offset) const;
  // Return the seqset entry for this specific readmap entry
  seqset_range read_ctx_entry(const readmap& rm, uint32_t readentry) const;
  // Makes an invalid end seqset range
  inline seqset_range end() const;
  // Find a sequence
  seqset_range find(const dna_slice &seq) const;
  seqset_range find(const dna_sequence &seq) const;
  // Finds an existing entry in the seqset given the precondition that
  // it exists.  For existing sequences, find_existing(seq) is
  // equivalent to find(seq).begin().  For nonexistant sequences, the
  // behavior of find_existing is undefined.
  uint64_t find_existing(const dna_slice& seq) const;
  // Same as find_existing, but will be faster if it's unique in the
  // first expected_unique_len bases.
  uint64_t find_existing_unique(const dna_slice& seq, size_t expected_unique_len) const;

  // Find inexact matches, returns false if more than max_results matches exist
  bool find_near(std::vector<seqset_range> &out, const dna_slice &seq,
                 size_t max_mismatch, size_t max_results) const;
  //bool find_diff(std::vector<seqset_range> &out, const dna_slice &seq,
 //                size_t max_diff, size_t max_results) const;

  bool entry_has_front(uint64_t entry, dna_base b) const { return m_prev[b]->get(entry); }
  uint64_t entry_push_front(uint64_t entry, dna_base b) const {
    return get_fixed(int(b)) + m_prev[b]->count(entry);
  }
  unsigned entry_size(uint64_t entry) const {
    return m_entry_sizes->get(entry);
  }
  unsigned entry_shared(uint64_t entry) const { return m_shared->get(entry); }
  // Returns the first base for an entry
  dna_base entry_get_base(uint64_t entry) const;
  uint64_t entry_pop_front(uint64_t entry) const;

  // Equivalent to find(seq).begin()
  boost::optional<uint64_t> find_kmer(const dna_slice& seq) const;

  // Equivalent to find(base + ctx_entry(seqset_id).sequence()
  //                                     .subseq(0, kmer_size - 1))
  boost::optional<uint64_t> kmer_push_front(uint64_t seqset_id,
                                            unsigned kmer_size,
                                            dna_base baes) const;

  void populate_pop_front_cache(
      progress_handler_t progress = null_progress_handler) const;
  void clear_pop_front_cache() const {
    delete[] m_pop_front_cache;
    m_pop_front_cache = NULL;
  }
  bool is_pop_front_cached() const { return m_pop_front_cache; }

  std::string uuid() const { return m_uuid; }

  // Compute summary table for entry_shared to speed up push_front_drop.
  void init_shared_lt_search();

private:
  struct open_gbwt_t {};
  void initialize_from_spiral_file(const spiral_file_open_state& state);

  struct seqset_metadata {
    TRANSFER_OBJECT {
      VERSION(0);
      FIELD(num_entries);
    };
    uint64_t num_entries = 0;
  };

  // Base can't be a dna_base because it can point one past the end of
  // the fixed array.
  void set_fixed(int base, uint64_t new_fixed) {
    CHECK_GE(base, 0);
    CHECK_LE(base, 4);
    reinterpret_cast<uint64_t *>(m_mutable_fixed.mutable_data())[base] =
        new_fixed;
  }
  uint64_t get_fixed(int base) const {
    DCHECK_GE(base, 0);
    DCHECK_LE(base, 4);
    return reinterpret_cast<const uint64_t *>(m_fixed.data())[base];
  }
  // Compute read_len
  void compute_read_len();

  std::string m_path;
  size_t m_entries = 0; // Number of 'seqs' (non-persistant)
  mutable_membuf m_mutable_fixed;
  membuf m_fixed;
  dna_base_array<std::unique_ptr<bitcount>> m_prev;
  std::unique_ptr<mutable_packed_varbit_vector> m_mutable_entry_sizes;
  std::unique_ptr<int_map_interface> m_entry_sizes;
  std::unique_ptr<mutable_packed_varbit_vector> m_mutable_shared;
  std::unique_ptr<int_map_interface> m_shared;
  std::unique_ptr<less_than_search> m_shared_lt_search;
  std::string m_uuid;

  bool m_is_final =
      false; // Set once the seqset is fully constructed and immutable.
  unsigned m_read_len = 0; // Computed read length

  // Value of element i of the vector is the entry that results from calling
  // pop_front on entry i.  Since we're storing 8 bytes per entry, the memory
  // usage here is potentially very large, so if the vector is empty,
  // pop_front
  // works by the slower (lg entry_count) algorithm.
  //
  // 5 uint8_t's are used for each element.
  // The first uint8_t is (value >> 32).
  // The next 4 uint8_ts are the low 32 bits of the value, in host order.
  mutable uint8_t *m_pop_front_cache = nullptr;

  // Membuf owning memory used to store old "gwt" format data.
  membuf m_old_format_buffer;
};

// Alias for "seqset".  Deprecated.
using seqset_file = seqset;

typedef std::unordered_map<uint64_t, unsigned> overlaps_t;

struct overlap_result_t {
  uint64_t seqset_id;
  unsigned overlap_bases;
};

// A seqset_range class holds a pair of iterators, m_begin and m_end, into its
// associated seqset and corresponds to a unique sequence which is the prefix of
// all the entries between m_begin and m_end.  The default seqset_range
// corresponds
// to the empty sequence and m_begin and m_end span the entire seqset.
class seqset_range : boost::totally_ordered<seqset_range> {
  friend class seqset;

public:
  bool operator<(const seqset_range &rhs) const {
    if (m_begin != rhs.m_begin)
      return m_begin < rhs.m_begin;
    if (m_end != rhs.m_end)
      // If the beginning is the same and the end is farther, the first is a prefix of the second.
      return m_end > rhs.m_end;
    return m_seq_size < rhs.m_seq_size;
  }
  bool operator==(const seqset_range &rhs) const {
    return m_begin == rhs.m_begin && m_end == rhs.m_end &&
           m_seq_size == rhs.m_seq_size;
  }
  // Get the number of bases of the range sequence
  unsigned size() const { return m_seq_size; }
  // Give the begining offset of this range
  uint64_t begin() const { return m_begin; }
  // Get the ending of this range
  uint64_t end() const { return m_end; }
  // Check if the range is valid
  bool valid() const { return m_begin < m_end; }
  // Returns how many bases we share with the previous context
  // of the same size, must be valid
  unsigned shared() const { return m_seqset->entry_shared(m_begin); }
  // Returns the number of bases shared with the given range
  unsigned shared_prefix_length(const seqset_range& rhs) const;
  // Returns the 'next' range of the same size, must be valid
  seqset_range next() const;
  // Add a dna base to the front, must be valid before this call
  seqset_range push_front(const dna_base &base) const;
  // Add a dna base to the front, dropping as much as needed to make it happen
  seqset_range push_front_drop(const dna_base &base, unsigned min_ctx = 0) const;
  // Pops a base to the front, must be valid before this call
  seqset_range pop_front() const;
  // Retuns a range expanded by dropping the final 'count' bases
  seqset_range pop_back(size_t count = 1) const;
  // Returns a range expanded by dropping down to 'new_size' bases.
  // Afterwards, .size() will be <= new_size.
  seqset_range truncate(size_t new_size) const;

  // Return the first base in the sequence, must be valid with non-zero size
  dna_base front() const;
  // Returns the whole dna_sequence associated with this range
  // or the first s bases if size is specified (error if size > m_entry_size)
  // This take log(n)*s where n is entries in pwbt table, and s is the size
  dna_sequence sequence(int size = -1) const;
  // Find all the reads (entries that are maximal) that end with
  // this range.  Find up to max_reads reads.  Return false if more reads
  // than max_reads were found, otherwise return true.
  bool find_maximal_prefix_reads(
      std::set<seqset_range> &results, uint32_t max_reads,
      unsigned min_overlap,
      const seqset_bitmap_base &read_bitmap = seqset_bitmap_true()) const;
  // Find all the reads (entries with length equal to read_size FROM READMAP)
  // that end with this range.  Find up to max_reads reads.  Return false
  // if more reads than max_reads were found, otherwise return true.
  bool find_full_prefix_reads(
      std::vector<seqset_range> &results, uint32_t max_reads,
      unsigned min_overlap,
      const readmap &read_bitmap) const;

  // Find all the reads (entries with length equal to read_size) that end with
  // this range sequence or a prefix of this range sequence at least min_overlap
  // in size.
  // Find up to max_reads reads.  If multiple output reads would overlap each
  // other,
  // return only the 'best' overlap, IE: only the most overlapping read for each
  // potential assembly is returned.  Return false if more reads than max_reads
  // ere found, otherwise return true.
  bool find_overlap_reads(
      overlaps_t &results, uint32_t max_reads, unsigned min_overlap,
      const seqset_bitmap_base &read_bitmap = seqset_bitmap_true(),
      bool rely_on_read_bitmap = false, unsigned added = 0) const;

  std::vector<overlap_result_t> find_overlap_reads_fair(
      uint32_t max_reads, unsigned min_overlap,
      const seqset_bitmap_base &read_bitmap = seqset_bitmap_true(),
      bool rely_on_read_bitmap = false, unsigned added = 0) const;

  // Returns true if this range is maximal, e.g. if there is no call
  // to push_front that returns a valid seqset_range and if size()
  // matches the size of the single entry this range refers to.
  bool is_maximal() const;
  // Returns true if there is a read in the readmap equal to the length
  bool is_full_read(const readmap &read_bitmap) const;
  // Returns true if this refers to a single full seqset entry.  This
  // check is less stringent than is_maximal; we may still be able to
  // push_front.
  bool is_seqset_entry() const {
    return begin() + 1 == end() && size() == m_seqset->entry_size(begin());
  }
  // Returns the single full seqset entry id associated with this range.
  uint64_t seqset_id() const {
    CHECK(is_seqset_entry());
    return begin();
  }

  // Bump the end pointer to the end of the seqset.
  void bump_to_end() { m_end = m_seqset->size(); }

  seqset_range() : m_seqset(nullptr), m_seq_size(0), m_begin(0), m_end(0) {}

  // Construct a range spanning the entire SEQSET
  seqset_range(const seqset *the_seqset_ptr)
      : m_seqset(the_seqset_ptr), m_seq_size(0), m_begin(0),
        m_end(the_seqset_ptr->size()) {}

  const seqset* get_seqset() const { return m_seqset; }

private:
  // Private constructor
  seqset_range(const seqset *seqset_, unsigned context_, uint64_t begin_,
               uint64_t end_)
      : m_seqset(seqset_), m_seq_size(context_), m_begin(begin_), m_end(end_) {}

  // Internal function that does work shared by seqset.begin() and next()
  void inner_next();
  // Internal version of pop that doesn't fix range length
  uint64_t inner_pop_front(dna_base b, uint64_t offset) const;

  const seqset *m_seqset;
  unsigned m_seq_size;
  uint64_t m_begin;
  uint64_t m_end;
};

inline seqset_range seqset::ctx_entry(uint64_t offset) const {
  if (offset >= size()) {
    throw(io_exception("Offset " + std::to_string(offset) + " larger than seqset size " +
                       std::to_string(size())));
  }
  return seqset_range(this, entry_size(offset), offset, offset + 1);
}

inline seqset_range seqset::end() const { return seqset_range(this, 0, 0, 0); }

inline uint64_t seqset::entry_pop_front(uint64_t entry) const {
  if (is_pop_front_cached()) {
    const uint8_t* cache_entry = &m_pop_front_cache[entry * 5];
    size_t hi = *cache_entry;
    hi <<= 32;
    ++cache_entry;
    return hi + *reinterpret_cast<const uint32_t*>(cache_entry);
  }

  return ctx_entry(entry).pop_front().begin();
}

struct seqset_range_hash {
  size_t operator()(const seqset_range& r) const {
    size_t val = r.begin();
    boost::hash_combine(val, r.end());
    return val;
  }
};

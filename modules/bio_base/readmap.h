#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_bitmap.h"
#include "modules/io/file_io.h"
#include "modules/io/int_map_interface.h"
#include "modules/io/packed_varbit_vector.h"
#include "modules/io/packed_vector.h"
#include "modules/io/sparse_multi.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/io/version.h"

#include <boost/range.hpp>

extern const product_version k_readmap_version;

struct alignas(8) readmap_header {
  uint64_t m_magic = k_magic;
  size_t m_offsets_offset = std::numeric_limits<size_t>::max();

  static constexpr uint64_t k_magic = 0x05D705D905E005DF;
};

// For old style readmap with gross/fine ids
struct readmap_info {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(version, TF_STRICT);
    FIELD(seqset_entry_count, TF_STRICT);
    FIELD(seqset_uuid, TF_STRICT);
    FIELD(user_string, TF_STRICT);
  }

  product_version version;
  uint64_t seqset_entry_count = 0;
  std::string seqset_uuid;
  std::string user_string;

  readmap_info() = default;
  explicit readmap_info(uint64_t the_seqset_entry_count)
      : version(k_readmap_version), seqset_entry_count(the_seqset_entry_count) {}
};

struct readmap_offsets {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(gross_read_ids, TF_STRICT);
    FIELD(fine_read_ids, TF_STRICT);
    FIELD(read_lengths, TF_STRICT);
    FIELD(info, TF_STRICT);
  }

  size_t gross_read_ids = 0;
  size_t fine_read_ids = 0;
  size_t read_lengths = 0;
  size_t info = 0;

  size_t gross_read_count() const {
    return (fine_read_ids - gross_read_ids) / sizeof(uint32_t) - 1;
  }
};

// For spiral file based readmap:
struct readmap_metadata {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(seqset_uuid);
  }

  std::string seqset_uuid;
};

class seqset;
class readmap : public seqset_bitmap_base {
  friend class make_readmap;

 public:
  struct pair_stats {
    size_t paired_reads = 0;
    size_t paired_bases = 0;
    size_t unpaired_reads = 0;
    size_t unpaired_bases = 0;
  };

  // Load a file created by make_readmap.  Caller is responsible for
  // making sure the_seqset continues to exist until the
  // readmap is destroyed.
  readmap(const std::shared_ptr<seqset>& the_seqset, const std::string& readmap_file_path,
          const spiral_file_options& sfopts);
  readmap(const std::shared_ptr<seqset>& the_seqset, const std::string& readmap_file_path)
      : readmap(the_seqset, readmap_file_path, m_spiral_file_opts) {}

  const std::string& path() const { return m_path; }

  // Open an anonymous readmap that's not associated with a seqset.
  // Deprecated.
  static std::unique_ptr<readmap> open_anonymous_readmap(const std::string& readmap_file_path);

  // API calls
  bool get_bit(uint64_t loc) const override;
  std::vector<int> fake_coverage(const dna_slice& seq);

  // "Approx" coverage does not properly calculate coverage in some
  // cases where reads are shorter than the length of seqset entries.
  //
  // Total coverage of sequence (fwd and rev)
  std::vector<int> approx_coverage(const dna_slice& seq) const;
  // If forward: returns only fwd coverage, else rev coverage
  std::vector<int> approx_strand_coverage(const dna_slice& seq, const bool& forward) const;
  // ret[0] is fwd coverage, ret[1] is rev coverage
  std::vector<std::vector<int>> approx_strand_coverage_split(const dna_slice& seq) const;

  std::pair<uint64_t, uint64_t> entry_to_index(uint64_t entry_id) const {
    return m_sparse_multi->lookup(entry_id);
  }
  uint64_t index_to_entry(uint64_t read_id) const {
    return m_sparse_multi->reverse_lookup(read_id);
  }
  std::pair<uint64_t, uint64_t> entry_to_index_range(uint64_t entry_id_start,
                                                     uint64_t entry_id_limit) const {
    return m_sparse_multi->lookup_range(entry_id_start, entry_id_limit);
  }

  size_t size() const { return m_read_lengths->size(); }

  // Support for pairing data:
  bool has_pairing_data() const { return m_pairing_data_present; }
  int get_readlength(uint32_t index) const;
  bool has_mate(uint32_t index) const;
  uint32_t get_mate(uint32_t index) const;
  uint64_t get_mate_entry(uint32_t index) const;
  bool get_is_forward(uint32_t index) const;
  // Returns the read id of the reverse complement of the given read.
  uint32_t get_rev_comp(uint32_t index) const;
  // Returns the reverse complement of the mate.  Mate loops must be
  // enabled.
  uint32_t get_mate_rc(uint32_t id) const;
  size_t get_read_count() const { return m_read_lengths->size() / 2; }
  size_t get_num_bases() const;
  pair_stats get_pair_stats() const;
  readmap_metadata metadata() const { return m_metadata; }

  // Enable mate loops; this allows for faster lookup for get_rev_comp.
  bool has_mate_loop() const { return m_mate_loop_ptr ? true : false; }
  void enable_mate_loop(
      std::function<dna_sequence(uint64_t /* seqset_id */, unsigned)> lookup_seq = {},
      progress_handler_t progress = null_progress_handler);
  unsigned min_read_len() const {
    calc_read_len_limits_if_needed();
    return m_min_read_len;
  }
  unsigned max_read_len() const {
    calc_read_len_limits_if_needed();
    return m_max_read_len;
  }

  static constexpr uint32_t k_null_index = std::numeric_limits<uint32_t>::max();

  class read;
  read get_read_by_id(uint32_t read_id) const;

  // Takes the result from get_mid_id() and returns a seqset id
  uint64_t mid_to_entry(uint64_t mid_id) const;

  // Iterating through prefix reads
  class read_iterator;
  using read_iterator_range = boost::iterator_range<read_iterator>;
  read_iterator_range get_prefix_reads(const seqset_range& r, int min_read_len = 0) const;
  boost::optional<uint32_t> get_longest_prefix_read_id(const seqset_range& r) const;
  boost::optional<read> get_longest_prefix_read(const seqset_range& r) const;

  // Iterating through containing reads
  class containing_read_iterator;
  using containing_read_iterator_range = boost::iterator_range<containing_read_iterator>;
  containing_read_iterator_range get_reads_containing(const seqset_range& r) const;

  const std::shared_ptr<seqset>& get_seqset() const { return m_seqset; }

  void calc_read_len_limits_if_needed() const {
    if (m_read_lengths_calculated.load(std::memory_order_relaxed)) {
      // Apparently std::call_once is not as fast as we'd prefer in some cases.
      return;
    }
    std::call_once(m_calc_read_len_limits_once,
                   [this]() { const_cast<readmap*>(this)->calc_read_len_limits(); });
  }

  // Returns a list of membufs to cache if memory caching is requested.
  membuf_cachelist membufs() const;

private:
  spiral_file_options m_spiral_file_opts;
  readmap(const std::string& readmap_file_path);

  void calc_read_len_limits();

  std::string m_path;
  std::shared_ptr<seqset> m_seqset;
  std::unique_ptr<sparse_multi> m_sparse_multi;

  // For new format:
  void open_spiral_file(const spiral_file_open_state& state);
  std::unique_ptr<spiral_file_open_mmap> m_opened;
  readmap_metadata m_metadata;

  std::atomic<bool> m_read_lengths_calculated {false};
  std::unique_ptr<int_map_interface> m_read_lengths;
  bool m_pairing_data_present = false;

  mutable std::once_flag m_calc_read_len_limits_once;
  unsigned m_min_read_len = std::numeric_limits<unsigned>::max();
  unsigned m_max_read_len = 0;

  std::unique_ptr<packed_vector<uint32_t, 32>> m_mate_pair_ptr;
  // mate loop is a pointer that runs in the following cyclical order for paired
  // reads:
  //
  //   Forward (is_forward = 1)
  //   RC (is_forward = 0)
  //   Pair (is_forward = 1)
  //   RC of pair (is_forward = 0)
  //
  // So, to get the reverse complement of a read, we follow this chain 1 time
  // if is_forward = 1, or 3 times if is_forward = 0.
  //
  // To get the pair, follow this chain twice.
  //
  // For non-paired reads, it runs in the following cyclical order:
  //
  //   Forward (is_forward = 1)
  //   RC (is_forward = 0)
  //
  // We can tell this is an unpaired read by following the chain twice
  // and returning to the original read.
  std::unique_ptr<int_map_interface> m_mate_loop_ptr;
  std::unique_ptr<packed_vector<unsigned, 1>> m_is_forward;
};

class readmap::read {
 public:
  read(const read&) = default;

  read get_mate() const { return read(m_readmap, m_readmap->get_mate(m_read_id)); }
  read get_rev_comp() const { return read(m_readmap, m_readmap->get_rev_comp(m_read_id)); }
  read get_mate_rc() const { return read(m_readmap, m_readmap->get_mate_rc(m_read_id)); }
  bool is_original_orientation() const { return m_readmap->get_is_forward(m_read_id); }
  int size() const { return m_readmap->get_readlength(m_read_id); }
  seqset_range get_seqset_entry() const {
    return m_readmap->m_seqset->ctx_entry(get_seqset_id()).truncate(size());
  }
  uint32_t get_read_id() const { return m_read_id; }
  uint64_t get_seqset_id() const {
    if (m_seqset_id == std::numeric_limits<uint64_t>::max()) {
      // Look up seqset id for this read if we haven't already.
      const_cast<read*>(this)->m_seqset_id = m_readmap->index_to_entry(m_read_id);
    }
    return m_seqset_id;
  }
  bool has_mate() const { return m_readmap->has_mate(m_read_id); }

  bool operator==(const read& rhs) const {
    return m_read_id == rhs.m_read_id && m_readmap == rhs.m_readmap;
  }
  bool operator!=(const read& rhs) const { return !(*this == rhs); }

  friend std::ostream& operator<<(std::ostream& os, const read& r);

  // Returns an id that is the same as another read's id if and only
  // if the seqset ids are the same.
  uint64_t get_mid_id() const {
    return m_readmap->m_sparse_multi->lookup_dest_to_mid(m_read_id);
  }

 private:
  read() = default;
  read(const readmap* rm)
      : m_readmap(rm) {}
  read(const readmap* rm, uint32_t read_id, uint64_t seqset_id)
      : m_readmap(rm), m_read_id(read_id), m_seqset_id(seqset_id) {}
  read(const readmap* rm, uint32_t read_id)
      : m_readmap(rm), m_read_id(read_id) {}
  friend class readmap;

  const readmap* m_readmap = nullptr;
  uint32_t m_read_id = std::numeric_limits<uint32_t>::max();
  mutable uint64_t m_seqset_id = std::numeric_limits<uint64_t>::max();
};

class readmap::read_iterator
    : public std::iterator<std::forward_iterator_tag, const readmap::read> {
 public:
  // Default is the end iterator.
  read_iterator() = default;
  read_iterator(const readmap* rm, uint32_t read_id, uint64_t seqset_id, int min_size, int max_size);

  const read& operator*() const
#ifdef __clang__
      // TODO(nils): Why does python readmap_test break without this when using clang?
      __attribute__((noinline))
#endif
  {
    return m_read;
  }
  read_iterator& operator++() {
    advance();
    skip_non_matching();
    return *this;
  }
  read_iterator operator++(int) {
    read_iterator result = *this;
    advance();
    skip_non_matching();
    return result;
  }
  bool operator==(const read_iterator& rhs) const {
    if (m_phase == phase::DONE) {
      return rhs.m_phase == phase::DONE;
    }
    if (m_phase != rhs.m_phase) {
      return false;
    }
    return m_read.m_read_id == rhs.m_read.m_read_id;
  }
  bool operator!=(const read_iterator& rhs) const { return !(*this == rhs); }

  friend std::ostream& operator<<(std::ostream& os, const read_iterator& it);

 private:
  friend class readmap;
  enum class phase { FORWARD, BACKWARD, DONE };
  friend std::ostream& operator<<(std::ostream& os, phase p);

  void advance();

  void skip_non_matching();
  // Returns true if we found an entry and thus don't have to search anymore.
  bool skip_non_matching_once();

  void done_direction();

  phase m_phase = phase::DONE;

  read m_read;
  int m_min_size = std::numeric_limits<int>::min();
  int m_max_read_len = std::numeric_limits<int>::max();

  // Values to reset to once we change direction.
  uint32_t m_orig_read_id = std::numeric_limits<uint32_t>::max();
  uint64_t m_orig_seqset_id = std::numeric_limits<uint64_t>::max();
  int m_orig_max_read_len = std::numeric_limits<int>::max();
};

class readmap::containing_read_iterator
    : public std::iterator<std::forward_iterator_tag, const std::pair<int, readmap::read>> {
 public:
  // Default is the end iterator.
  containing_read_iterator() = default;
  containing_read_iterator(const readmap* rm, const seqset_range& r);

  const std::pair<int, readmap::read>& operator*() const { return m_offset_and_read; }
  containing_read_iterator& operator++() {
    CHECK(!at_end());
    advance_read();
    skip_non_matching();
    return *this;
  }
  containing_read_iterator operator++(int) {
    CHECK(!at_end());
    containing_read_iterator result = *this;
    advance_read();
    skip_non_matching();
    return result;
  }
  bool operator==(const containing_read_iterator& rhs) const {
    if (m_range != rhs.m_range) {
      return false;
    }
    if (!m_range.valid()) {
      return true;
    }
    return m_offset_and_read == rhs.m_offset_and_read;
  }
  bool operator!=(const containing_read_iterator& rhs) const { return !(*this == rhs); }

  friend std::ostream& operator<<(std::ostream& os, const containing_read_iterator& it);

 private:
  void advance_read();
  void skip_non_matching();
  void start_entry();
  void advance_entry();
  bool at_end() const;
  read& get_read() { return m_offset_and_read.second; }
  const read& get_read() const { return m_offset_and_read.second; }

  std::pair<int, read> m_offset_and_read{0, read()};
  seqset_range m_range;
  unsigned m_orig_len = 0;

  // Current range
  uint32_t m_end_read_id = std::numeric_limits<uint32_t>::max();
};

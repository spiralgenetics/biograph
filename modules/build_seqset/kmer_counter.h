#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"
#include "modules/build_seqset/kmer_count_table.h"
#include "modules/io/packed_vector.h"
#include "modules/io/string_view.h"
#include "modules/io/track_mem.h"

#include <boost/optional.hpp>

class kmer_set;

namespace build_seqset {

struct count_kmer_options {
  // Number of kmer partitions.  Writes to each of these partitions is
  // batched, so partitions should be small enough that random access
  // to it doesn't thrash the TLB cache.
  unsigned partitions = 256;

  // Size of kmers to use, in bases.  This must be less than 32.
  unsigned kmer_size = 30;

  // Minimum number of times a kmer must occur to prevent being filtered.
  unsigned min_count = 3;

  // Calculate number of passes so that we don't exceed this much
  // memory use per pass.
  size_t max_memory_bytes = 20ULL * 1024 * 1024;

  // If non-zero, the maximum number of entries that should be in the
  // probabilistic table.  This makes it so we can avoid huge memory
  // allocations when processing small datasets.
  size_t max_prob_table_entries = 0;

  // Don't exceed this density for exact count tables.
  double max_exact_table_density = 0.7;

  // Don't go lower than this density for exact count tables.
  double min_exact_table_density = 0.1;

  // Portion of (post-filtered) kmers that we expect to see more than
  // 255 times.
  double overflow_table_size_ratio = 0.05;

  // Absolute minimum number of entries in a table, even in small
  // datasets.
  size_t abs_min_table_size = 512ULL * 1024;

  // Number of entries in the probabilistic kmer table.  Each of these
  // entries takes 2 bits of RAM during the probablistic phase, and 1
  // bit of RAM during the exact count phase.  If 0, this is
  // calculated based on available memory.
  size_t prob_table_entries = 0;

  // Number of operations to be queued up before writing to a
  // partition.
  unsigned partition_batch_size = 256;

  // Keep temporary files instead of deleting them when done.
  bool keep_temporaries = false;

  // Force a certain number of exact passes.  Zero means autodetec
  // based on memory available.
  unsigned force_exact_passes = 0;

  // Progress handler
  progress_handler_t progress = null_progress_handler;

  static const count_kmer_options& defaults();
};

// kmer_counter implements a 2-stage count of "kmer_t"s, each stage
// which has 1 or more passes over the data.
//
// On the "prob" passes, a probabilistic table of 2 bit counters is
// filled based on hashes of the kmers.
//
// On the "exact" passes, an exact table of counters is counted for
// each kmer.
//
// To use kmer_counter, follow this control flow:
//
// kmer_counter counter(options);
// counter.start_prob_pass();
//   In parallel, execute:
//     kmer_counter::prob_pass_processor p(counter);
//     p.add(sequence1)
//     p.add(sequence2)
//     ...
// counter.close_prob_pass();
// For pass in [0, exact_passes):
//   counter.start_exact_pass(pass);
//   In parallel, execute:
//     kmer_counter::exact_pass_processor p(counter);
//     p.add(sequence1)
//     p.add(sequence2)
//     ...
// counter.close_exact_passes();
// counter.extract_exact_counts(
//     [](kmer_counter::iterator start, kmer_counter::iterator limit) {
//        for (auto it = start; it != limit; ++it) {
//           if (it->fwd_count > threshold || it->rev_cout > threshold) {
//               printf("%ld passed the threshold!", it->kmer);
//           }
//     }, progress);
// counter.close();
//
// You may use multiple pass_processors in parallel for the same
// kmer_counter, but you may not call .add() in parallel on the same
// pass_processor.
//
// extract_exact_counts may call the supplied std::function in parallel,
// so it must be thread safe.
//
// You may call extract_exact_counts multiple times; each time
// it's called it will iterate through all the kmers.
//
// extract_exact_counts is not guaranteed to return any entries with
// counts less than 3, nor is it guaranteed not to.

class kmer_counter {
 public:
  kmer_counter(
      const count_kmer_options& options = count_kmer_options::defaults());

  struct element {
    kmer_t kmer = std::numeric_limits<kmer_t>::max();
    uint32_t fwd_count = 0;
    uint32_t rev_count = 0;

    // True if any reads start with this kmer.
    bool fwd_starts_read = false;

    // True if any reads end with this kmer.
    bool rev_starts_read = false;
  };

  class extract_iterator;
  class prob_pass_processor;
  class exact_pass_processor;

  unsigned exact_passes() const { return m_exact_passes; }
  void start_prob_pass();
  void close_prob_pass();

  void start_exact_pass(unsigned pass_num);
  void close_exact_passes();
  void extract_exact_counts(
      const std::function<void(extract_iterator, extract_iterator)>& output_f
  );
  void close();

  void set_progress_handler(progress_handler_t progress);

 private:
  class pass_processor;

  void create_prob_filters();
  void create_exact_counters(unsigned pass_num);

  void close_exact_pass();
  void show_exact_stats();

  // Returns which partition is associated with the given kmer.
  unsigned kmer_partition(kmer_t kmer) const {
    // We have to add mixing steps here instead of just multiplying;
    // otherwise it won't be evenly distributed if we don't have a
    // prime number of partitions.
    uint64_t hash = kmer;
    hash *= 0xff51afd7ed558ccd;
    hash ^= (hash >> 33);
    hash *= 0xc4ceb9fe1a85ec53;
    hash ^= (hash >> 33);
    uint64_t quotient = m_partition_divider.perform_divide(hash);
    return hash - quotient * m_options.partitions;
  }

  // Returns a hashed value of kmer.
  static uint64_t pt_hash_kmer(kmer_t kmer) {
    return kmer * 11304120250909662091ULL;
  }

  bool partition_is_active(unsigned pass_num, unsigned tot_passes,
                           unsigned part_num) const;

  // Tweaks approx_size so that it shares no factors with m_options.partitions.
  size_t get_relatively_prime_partition_size(size_t approx_size) const;

  count_kmer_options m_options;

  enum class count_state {
    INITIALIZED,
    PROB_PASS,
    PROB_PASS_FINISHED,
    EXACT_PASS,
    EXACT_PASSES_FINISHED,
    CLOSED
  };
  friend std::ostream& operator<<(std::ostream&, count_state);

  count_state m_count_state = count_state::INITIALIZED;
  unsigned m_pass_num = 0;
  unsigned m_exact_passes = 0;
  size_t m_exact_entries_needed = 0;

  std::string m_temp_dir;

  // Bitmask for kmers to trim them down to kmer length.
  kmer_t m_kmer_mask;

  // Probabilistic kmer count tables, per partition.
  using prob_table_t = packed_vector<unsigned, 1>;
  std::vector<boost::optional<prob_table_t>> m_prob_table;
  using mutable_prob_table_t = mutable_packed_vector<unsigned, 2>;
  std::vector<boost::optional<mutable_prob_table_t>> m_mutable_prob_table;

  // Exact kmer count tables, per partition.
  using exact_count_table_t = kmer_count_table<uint8_t>;
  std::vector<boost::optional<exact_count_table_t>> m_exact_table;

  // Overflow counts for all partitions.
  using exact_overflow_count_table_t = kmer_count_table<uint32_t>;
  boost::optional<exact_overflow_count_table_t> m_exact_overflow_table;

  // Divides by m_option.partitions
  libdivide::divider<uint64_t> m_partition_divider;

  std::atomic<size_t> m_prob_skipped{0};
  std::atomic<size_t> m_tot_exact_kmers{0};
};

class kmer_counter::extract_iterator
    : public boost::iterator_facade<extract_iterator, element const,
                                    std::random_access_iterator_tag, element> {
 public:
  extract_iterator(exact_count_table_t::const_iterator pos,
                   const exact_overflow_count_table_t* overflow_table)
      : m_pos(pos), m_overflow_table(overflow_table) {}

  element dereference() const {
    element result;
    result.kmer = m_pos->kmer();
    result.fwd_count = m_pos->fwd_count;
    result.rev_count = m_pos->rev_count;
    result.fwd_starts_read = m_pos->fwd_flag();
    result.rev_starts_read = m_pos->rev_flag();
    const auto& of = m_overflow_table->get(m_pos->kmer());
    if (of) {
      result.fwd_count += of.fwd_count;
      result.rev_count += of.rev_count;
    }
    return result;
  }
  bool equal(const extract_iterator& rhs) const { return m_pos == rhs.m_pos; }
  void advance(int64_t diff) { m_pos += diff; }
  void increment() { ++m_pos; }
  void decrement() { --m_pos; }
  int64_t distance_to(const extract_iterator& rhs) const {
    return rhs.m_pos - m_pos;
  }

 private:
  exact_count_table_t::const_iterator m_pos;
  const exact_overflow_count_table_t* m_overflow_table = nullptr;
};

class kmer_counter::pass_processor {
 public:
  static constexpr kmer_t k_fwd_flag = 1ULL << 63;
  static constexpr kmer_t k_rev_flag = 1ULL << 62;
  static constexpr kmer_t k_kmer_mask = (1ULL << 62) - 1;

  void add_kmer(kmer_t kmer, bool fwd_flag, bool rev_flag) {
    DCHECK_LE(kmer, m_kmer_counter.m_kmer_mask);
    kmer_t canon = canonicalize(kmer, m_kmer_size);
    unsigned part_num = m_kmer_counter.kmer_partition(canon);
    auto& bounds = m_part_bounds[part_num];
    auto& part_cur = bounds.first;
    auto& part_end = bounds.second;
    if (!part_cur) {
      return;
    }
    if (part_cur == part_end) {
      kmer_t* start = m_part_queues[part_num].data();
      flush_part(part_num, start, part_cur);
      part_cur = start;
    }
    *part_cur =
        kmer | (fwd_flag ? k_fwd_flag : 0) | (rev_flag ? k_rev_flag : 0);
    ++part_cur;
  }

  static std::atomic<void *> g_annotate_value;

  void add(string_view seq) {
    auto it = seq.begin();
    auto end = seq.end();

    bool is_first_kmer = true;

    unsigned initial_kmer_left = m_kmer_size;
    kmer_t kmer = 0;
    while (it != end) {
      if (*it == 'N') {
        ++it;
        initial_kmer_left = m_kmer_size;
        is_first_kmer = false;
        continue;
      }
      kmer <<= 2;
      kmer |= int(dna_base(*it));
      ++it;
      if (initial_kmer_left) {
        --initial_kmer_left;
      }

      if (!initial_kmer_left) {
        bool is_last_kmer = (it == end);
        add_kmer(kmer & m_kmer_mask, is_first_kmer,
                 is_last_kmer);
        is_first_kmer = false;
      }
    }
  }

  void flush_all();

 protected:
  pass_processor(kmer_counter& k, unsigned tot_passes);
  pass_processor& operator=(const pass_processor&) = delete;
  pass_processor(const pass_processor&) = delete;
  virtual ~pass_processor();

  virtual void flush_part(unsigned part_num, const kmer_t* start,
                          const kmer_t* limit) = 0;

 protected:
  kmer_counter& m_kmer_counter;
  unsigned m_pass_num;
  unsigned m_tot_passes;
  tracked_vector<std::vector<kmer_t>> m_part_queues;
  std::vector<std::pair<kmer_t* /* current */, kmer_t* /* end */>> m_part_bounds;
  const kmer_t m_kmer_mask;
  const unsigned m_kmer_size;
};

class kmer_counter::prob_pass_processor final : public pass_processor {
 public:
  prob_pass_processor(kmer_counter& k);
  ~prob_pass_processor();

 private:
  void flush_part(unsigned part_num, const kmer_t* start,
                  const kmer_t* limit) override;
};

class kmer_counter::exact_pass_processor final : public pass_processor {
 public:
  exact_pass_processor(kmer_counter& k);
  ~exact_pass_processor();

 private:
  void flush_part(unsigned part_num, const kmer_t* start,
                  const kmer_t* limit) override;
};

}  // namespace build_seqset

#include "modules/build_seqset/kmer_counter.h"
#include "modules/io/config.h"
#include "modules/io/parallel.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/io/stats.h"
#include "modules/io/track_mem.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <numeric>

namespace build_seqset {

std::ostream& operator<<(std::ostream& os, kmer_counter::count_state state) {
  switch (state) {
    case kmer_counter::count_state::INITIALIZED:
      os << "INITIALIZED";
      return os;
    case kmer_counter::count_state::PROB_PASS:
      os << "PROB_PASS";
      return os;
    case kmer_counter::count_state::PROB_PASS_FINISHED:
      os << "PROB_PASS_FINISHED";
      return os;
    case kmer_counter::count_state::EXACT_PASS:
      os << "EXACT_PASS";
      return os;
    case kmer_counter::count_state::EXACT_PASSES_FINISHED:
      os << "EXACT_PASSES_FINISHED";
      return os;
    case kmer_counter::count_state::CLOSED:
      os << "CLOSED";
      return os;
  };
  os << "Unknown count state " << int(state) << "\n";
  return os;
}

const count_kmer_options& count_kmer_options::defaults() {
  static count_kmer_options g_default_options;
  return g_default_options;
}

kmer_counter::kmer_counter(const count_kmer_options& options)
    : m_options(options),
      m_kmer_mask(~(std::numeric_limits<kmer_t>::max() << (options.kmer_size * 2))),
      m_prob_table(m_options.partitions),
      m_mutable_prob_table(m_options.partitions),
      m_exact_table(m_options.partitions),
      m_partition_divider(m_options.partitions) {
  if (options.kmer_size > 31) {
    throw(io_exception("A maximum kmer size of 31 is supported for read correction"));
  }
  m_temp_dir = CONF_S(temp_root);
  if (boost::starts_with(m_temp_dir, "file://")) {
    m_temp_dir = m_temp_dir.substr(7);
  }
  boost::filesystem::create_directories(m_temp_dir);

  size_t prob_table_entries = m_options.max_memory_bytes *
                              4
                              // Leave room for a partition to be compressed (2 bits per entry
                              // to one bit per entry)
                              * m_options.partitions / (m_options.partitions + 1);

  if (m_options.max_prob_table_entries && prob_table_entries > m_options.max_prob_table_entries) {
    SPLOG("Limiting probabilistic table entries from %ld to %ld", prob_table_entries,
          m_options.max_prob_table_entries);
    prob_table_entries = m_options.max_prob_table_entries;
  }
  if (m_options.prob_table_entries) {
    SPLOG("Overriding probabilistic table entries from %ld to user-supplied %ld",
          m_options.prob_table_entries, prob_table_entries);
  } else {
    SPLOG("Using %ld probabilistic table entries, %.2f MB RAM", prob_table_entries,
          prob_table_entries / 4. / 1024. / 1024.);
    m_options.prob_table_entries = prob_table_entries;
  }

  if (m_options.prob_table_entries < m_options.abs_min_table_size) {
    SPLOG("Increasing probabilistic table size from %ld to absolute minimum %ld",
          m_options.prob_table_entries, m_options.abs_min_table_size);
    m_options.prob_table_entries = m_options.abs_min_table_size;
  }
}

void kmer_counter::create_prob_filters() {
  size_t prob_filter_size =
      get_relatively_prime_partition_size(m_options.prob_table_entries / m_options.partitions);
  SPLOG("Creating probabilistic filters with %ld entries per partition, %d partitions",
        prob_filter_size, m_options.partitions);

  parallel_for(
      0, m_options.partitions,
      [&](size_t part_num) {
        m_mutable_prob_table[part_num].emplace(prob_filter_size, "kmer_counter_build_prob");
      },
      m_options.progress);
}

void kmer_counter::create_exact_counters(unsigned pass_num) {
  CHECK_LT(pass_num, m_exact_passes);

  unsigned active_partitions = 0;
  for (unsigned part_num = 0; part_num < m_options.partitions; ++part_num) {
    if (partition_is_active(pass_num, m_exact_passes, part_num)) {
      active_partitions++;
    }
  }
  CHECK(active_partitions);

  size_t overflow_table_bytes =
      m_exact_overflow_table->size() * sizeof(exact_overflow_count_table_t::element);
  CHECK_LT(overflow_table_bytes, m_options.max_memory_bytes);
  size_t memory_per_partition =
      (m_options.max_memory_bytes - overflow_table_bytes) / active_partitions;
  size_t prob_memory_per_partition = m_options.prob_table_entries / m_options.partitions / 8;
  size_t exact_memory_per_partition = memory_per_partition - prob_memory_per_partition;

  double bytes_per_exact_entry = sizeof(exact_count_table_t::element);

  size_t max_exact_entries = m_exact_entries_needed * active_partitions *
                             m_options.max_exact_table_density / m_options.partitions /
                             m_options.min_exact_table_density;
  size_t exact_entries = exact_memory_per_partition / bytes_per_exact_entry;
  if (exact_entries > max_exact_entries) {
    SPLOG(
        "Small dataset; decreasing exact entries per partition from %ld to "
        "%ld.",
        exact_entries, max_exact_entries);
    exact_entries = max_exact_entries;
  }
  if (exact_entries < m_options.abs_min_table_size) {
    SPLOG("Increasing exact table size from %ld to absolute minimum %ld", exact_entries,
          m_options.abs_min_table_size);
    exact_entries = m_options.abs_min_table_size;
  }
  size_t exact_table_size = get_relatively_prime_partition_size(exact_entries);

  SPLOG(
      "Creating exact counter with %ld entries for %d "
      "partitions.  RAM use: %.2f MB for prob table, %.2f MB for exact table",
      exact_table_size, active_partitions,
      prob_memory_per_partition * active_partitions / 1024. / 1024,
      exact_table_size * sizeof(exact_count_table_t::element) * active_partitions / 1024. / 1024);

  parallel_for(  //
      0, m_options.partitions, [this, pass_num, exact_table_size](size_t part_num) {
        if (partition_is_active(pass_num, m_exact_passes, part_num)) {
          spiral_file_options opts;
          opts.read_into_ram = true;
          // Open the probabilistic table for filtering.
          spiral_file_open_mmap o(m_temp_dir + "/kmerize_prob-part-" + std::to_string(part_num));
          m_prob_table[part_num].emplace(o.open(), "kmer_counter_prob");
          m_exact_table[part_num].emplace(exact_table_size, "main");
        }
      });
}

namespace {
size_t gcd(size_t a, size_t b) {
  if (b == 0) {
    return a;
  }
  return gcd(b, a % b);
}

}  // namespace

size_t kmer_counter::get_relatively_prime_partition_size(size_t approx_size) const {
  size_t result = approx_size;
  while (gcd(result, m_options.partitions) > 1) {
    ++result;
  }
  return result;
}

bool kmer_counter::partition_is_active(unsigned pass_num, unsigned tot_passes,
                                       unsigned part_num) const {
  unsigned processed_in_pass = part_num * tot_passes / m_options.partitions;
  return pass_num == processed_in_pass;
}

kmer_counter::pass_processor::pass_processor(kmer_counter& k, unsigned tot_passes)
    : m_kmer_counter(k),
      m_pass_num(m_kmer_counter.m_pass_num),
      m_tot_passes(tot_passes),
      m_part_queues(k.m_options.partitions, track_alloc("kmer_counter:part_queue")),
      m_part_bounds(k.m_options.partitions),
      m_kmer_mask(k.m_kmer_mask),
      m_kmer_size(k.m_options.kmer_size) {
  CHECK_LT(m_pass_num, m_tot_passes);
  for (unsigned i = 0; i < k.m_options.partitions; ++i) {
    if (m_kmer_counter.partition_is_active(m_pass_num, tot_passes, i)) {
      m_part_queues[i].resize(k.m_options.partition_batch_size);
      m_part_bounds[i] = std::make_pair(m_part_queues[i].data(),
                                        m_part_queues[i].data() + m_part_queues[i].size());
    }
  }
}

void kmer_counter::pass_processor::flush_all() {
  for (unsigned i = 0; i < m_part_queues.size(); ++i) {
    auto& bounds = m_part_bounds[i];
    auto& part_cur = bounds.first;
    if (part_cur) {
      flush_part(i, m_part_queues[i].data(), part_cur);
      part_cur = m_part_queues[i].data();
    }
  }
}

kmer_counter::pass_processor::~pass_processor() {
  for (unsigned i = 0; i < m_part_queues.size(); ++i) {
    auto& bounds = m_part_bounds[i];
    auto& part_cur = bounds.first;
    if (part_cur) {
      CHECK(m_part_queues[i].data() == part_cur)
          << "Must call flush_all before destroying pass processor";
    }
  }
}

void kmer_counter::start_prob_pass() {
  CHECK_EQ(m_count_state, count_state::INITIALIZED);
  m_count_state = count_state::PROB_PASS;

  create_prob_filters();
  SPLOG("kmer_counter: starting probabilistic pass");
  track_mem::reset_stats();
}

void kmer_counter::close_prob_pass() {
  CHECK_EQ(m_count_state, count_state::PROB_PASS);

  SPLOG("Closing probalistic pass");
  constexpr unsigned k_histo_size = 1 + mutable_prob_table_t::max_value_static();
  static constexpr bool k_show_part_stats = false;
  std::array<simple_stats<double>, k_histo_size> part_histo{};
  std::array<size_t, k_histo_size> histo{};

  unsigned min_count = mutable_prob_table_t::max_value_static();
  if (m_options.min_count < min_count) {
    min_count = m_options.min_count;
    CHECK_GT(min_count, 0);
  }

  std::mutex mu;
  size_t tot_entries = 0;
  size_t prob_bytes = m_options.prob_table_entries / 4;
  size_t prob_output_bytes_per_partition =
      (m_options.prob_table_entries + (m_options.partitions * 8 - 1)) / (m_options.partitions * 8);

  size_t parts_in_parallel = 0;
  size_t free_mem = 0;
  if (m_options.max_memory_bytes > prob_bytes) {
    free_mem = m_options.max_memory_bytes - prob_bytes;
  }
  parts_in_parallel = free_mem / prob_output_bytes_per_partition;

  if (parts_in_parallel < 1) {
    SPLOG("Warning: May not have enough memory to save probablistic pass data; free_mem = %ld",
          free_mem);
    parts_in_parallel = 1;
  } else {
    SPLOG(
        "Saving probablistic table entries starting out %ld partitions at once since free mem = "
        "%ld and bytes per partition = %ld",
        parts_in_parallel, free_mem, prob_output_bytes_per_partition);
  }

  std::condition_variable cv;

  parallel_for(
      0, m_options.partitions,
      [&](size_t part_num) {
        if (!m_mutable_prob_table[part_num]) {
          return;
        }

        {
          std::unique_lock<std::mutex> l(mu);
          cv.wait(l, [&parts_in_parallel]() { return parts_in_parallel > 0; });
          --parts_in_parallel;
        }

        spiral_file_options sfopts;
        const auto& t = *m_mutable_prob_table[part_num];
        std::array<size_t, k_histo_size> phisto{};
        std::string prob_filename = m_temp_dir + "/kmerize_prob-part-" + std::to_string(part_num);
        {
          spiral_file_create_mmap c(prob_filename, sfopts.with_delayed_write(true));
          mutable_packed_vector<unsigned, 1> output(c.create(), t.size());

          for (size_t i = 0; i != t.size(); ++i) {
            auto val = t[i];
            phisto[val]++;
            if (val >= min_count) {
              output.at(i).set_unlocked(1);
            }
          }
        }

        std::lock_guard<std::mutex> l(mu);
        tot_entries += t.size();
        for (unsigned i = 0; i < k_histo_size; ++i) {
          part_histo[i].add_sample(phisto[i] * 100. / t.size());
          histo[i] += phisto[i];
        }
        m_mutable_prob_table[part_num].reset();
        parts_in_parallel +=
            // We freed up one input partition which is twice the size of an output partition:
            2
            // Plus one output partition:
            + 1;
        cv.notify_all();
        if (k_show_part_stats) {
          SPLOG(
              "Partition %ld count: total: %ld 0: %ld (%.2f%%) 1: %ld (%.2f%%) 2: "
              "%ld (%.2f%%) 3: %ld (%.2f%%) ",
              part_num, t.size(), phisto[0], phisto[0] * 100. / t.size(), phisto[1],
              phisto[1] * 100. / t.size(), phisto[2], phisto[2] * 100. / t.size(), phisto[3],
              phisto[3] * 100. / t.size());
        }
      },
      m_options.progress);

  SPLOG("%ld probabilistic bitmap entries with the following counts:", tot_entries);
  size_t kmers_set = 0;
  size_t passing_kmers_set = 0;
  for (unsigned i = 0; i < k_histo_size; ++i) {
    auto& h = part_histo[i];
    h.analyze();
    if (i > 0) {
      kmers_set += histo[i];
    }
    if (i >= min_count) {
      passing_kmers_set += histo[i];
    }
    SPLOG(" %d: %12ld (%6.2f%% avg per partition, %6.2f%% min, %6.2f%% max)", i, histo[i], h.avg,
          h.min, h.max);
  }

  double kmers_ratio = kmers_set * 1. / tot_entries;
  double passing_kmers_ratio = passing_kmers_set * 1. / tot_entries;
  CHECK_GE(kmers_ratio, passing_kmers_ratio);

  double est_kmers = kmers_set / (1 - kmers_ratio);
  double est_passing_kmers = passing_kmers_set / (1 - passing_kmers_ratio);
  CHECK_GE(est_kmers, est_passing_kmers);

  SPLOG("Estimating %ld total kmers, %ld (%.2f%%) of which need exact counts.", size_t(est_kmers),
        size_t(est_passing_kmers), est_passing_kmers * 100 / est_kmers);

  double exact_table_entries_used = est_passing_kmers +
                                    // False positives in the probabilistic filter:
                                    (est_kmers - est_passing_kmers) * passing_kmers_ratio;

  size_t overflow_table_size = exact_table_entries_used * m_options.overflow_table_size_ratio /
                               m_options.max_exact_table_density;
  if (overflow_table_size < m_options.abs_min_table_size) {
    SPLOG("Increasing overflow table size from %ld to absolute minimum %ld", overflow_table_size,
          m_options.abs_min_table_size);
    overflow_table_size = m_options.abs_min_table_size;
  }

  size_t overflow_table_bytes = overflow_table_size * sizeof(exact_overflow_count_table_t::element);
  SPLOG("Overflow table has %ld entries using %.2f MB RAM", overflow_table_size,
        overflow_table_bytes / 1024. / 1024);
  m_exact_overflow_table.emplace(overflow_table_size, "kmer_counter_overflow_table");

  m_exact_entries_needed = exact_table_entries_used / m_options.max_exact_table_density;
  SPLOG("Requiring %ld exact table entries", m_exact_entries_needed);

  if (overflow_table_bytes >= m_options.max_memory_bytes) {
    SPLOG(
        "Warning: overflow table using %ld bytes, which is more than max memory setting %ld; "
        "increasing max memory setting.",
        overflow_table_bytes, m_options.max_memory_bytes);
    m_options.max_memory_bytes = overflow_table_bytes + 1;
  }

  size_t total_exact_size =
      m_options.prob_table_entries / 8 +  // 8 prob table entries per byte after built
      m_exact_entries_needed * sizeof(exact_count_table_t::element);
  m_exact_passes = (total_exact_size / (m_options.max_memory_bytes - overflow_table_bytes)) + 1;
  if (m_exact_passes > m_options.partitions) {
    SPLOG("Limiting exact passes from %d to the number of partitions, %d", m_exact_passes,
          m_options.partitions);
    m_exact_passes = m_options.partitions;
  }

  if (m_options.force_exact_passes) {
    SPLOG("Overriding exact passes to %d from %d.", m_options.force_exact_passes, m_exact_passes);
    m_exact_passes = m_options.force_exact_passes;
  }

  SPLOG("Exact entries need %.2f MB of memory; using %d passes.", total_exact_size / 1024. / 1024,
        m_exact_passes);

  m_options.progress(1.0);
  m_count_state = count_state::PROB_PASS_FINISHED;
}

void kmer_counter::start_exact_pass(unsigned pass_num) {
  if (pass_num == 0) {
    CHECK_EQ(m_count_state, count_state::PROB_PASS_FINISHED);
  } else {
    CHECK_EQ(m_count_state, count_state::EXACT_PASS);
    CHECK_EQ(m_pass_num + 1, pass_num);
    close_exact_pass();
  }
  m_pass_num = pass_num;
  m_count_state = count_state::EXACT_PASS;
  m_prob_skipped.store(0);
  m_tot_exact_kmers.store(0);

  create_exact_counters(pass_num);
  SPLOG("kmer_counter: starting exact pass %d/%d.", pass_num + 1, m_exact_passes);
  track_mem::reset_stats();
}

void kmer_counter::close_exact_pass() {
  CHECK_EQ(m_count_state, count_state::EXACT_PASS);
  show_exact_stats();
  SPLOG("Saving exact counts");
  parallel_for(0, m_options.partitions, [&](size_t part_num) {
    if (!m_exact_table[part_num]) {
      return;
    }
    CHECK(m_prob_table[part_num]);
    m_prob_table[part_num].reset();

    if (!m_options.keep_temporaries) {
      std::string prob_file = m_temp_dir + "/kmerize_prob-part-" + std::to_string(part_num);
      boost::filesystem::remove(prob_file);
    }

    auto& et = m_exact_table[part_num];
    spiral_file_options sfopts;
    std::string exact_file = m_temp_dir + "/kmerize_exact-part-" + std::to_string(part_num);
    spiral_file_create_mmap c(exact_file, sfopts.with_delayed_write(true));
    et->save(c.create());
    et.reset();
  });
  SPLOG("Done saving exact counts");
}

void kmer_counter::close_exact_passes() {
  CHECK_EQ(m_count_state, count_state::EXACT_PASS);
  CHECK_EQ(m_pass_num + 1, m_exact_passes);
  close_exact_pass();

  for (unsigned part_num = 0; part_num < m_options.partitions; ++part_num) {
    if (!m_exact_table[part_num]) {
      std::string exact_file = m_temp_dir + "/kmerize_exact-part-" + std::to_string(part_num);
      spiral_file_open_mmap o(exact_file);
      m_exact_table[part_num].emplace(o.open());
      if (!m_options.keep_temporaries) {
        boost::filesystem::remove(exact_file);
      }
    }
  }

  size_t tot_overflow_used = 0;
  for (auto it = m_exact_overflow_table->begin(); it != m_exact_overflow_table->end(); ++it) {
    if (*it) {
      tot_overflow_used++;
    }
  }

  SPLOG("Overflow used: %ld/%ld (%.2f%%)", tot_overflow_used, m_exact_overflow_table->size(),
        tot_overflow_used * 100. / m_exact_overflow_table->size());

  m_count_state = count_state::EXACT_PASSES_FINISHED;
}

void kmer_counter::show_exact_stats() {
  CHECK_EQ(m_count_state, count_state::EXACT_PASS);
  SPLOG(
      "Exact pass %d/%d complete; %ld/%ld (%.2f%%) skipped due to "
      "probabilistic lookup",
      m_pass_num + 1, m_exact_passes, m_prob_skipped.load(), m_tot_exact_kmers.load(),
      m_prob_skipped.load() * 100. / m_tot_exact_kmers.load());

  std::mutex mu;
  size_t tot_exact_entries = 0;
  size_t tot_exact_entries_used = 0;
  simple_stats<double> part_used;

  constexpr unsigned k_histo_size = 3;
  std::array<std::atomic<size_t>, k_histo_size> tot_histo{};

  parallel_for(0, m_options.partitions, [&](size_t part_num) {
    if (!m_exact_table[part_num]) {
      return;
    }
    auto& et = *m_exact_table[part_num];

    size_t orig_size = et.size();
    et.compact();

    std::array<size_t, k_histo_size> local_histo{};

    for (auto it = et.begin(); it != et.end(); ++it) {
      size_t tot_count = it->fwd_count + it->rev_count;
      if (tot_count < local_histo.size()) {
        local_histo[tot_count]++;
      }
    }
    size_t part_exact_entries_used = et.end() - et.begin();

    for (unsigned i = 0; i < k_histo_size; ++i) {
      tot_histo[i].fetch_add(local_histo[i]);
    }

    {
      std::lock_guard<std::mutex> l(mu);
      part_used.add_sample(part_exact_entries_used * 100. / orig_size);
      tot_exact_entries += orig_size;
      tot_exact_entries_used += part_exact_entries_used;
    }
  });

  SPLOG("Exact pass main counters used %ld of %ld total entries (%.2f%%) ", tot_exact_entries_used,
        tot_exact_entries, tot_exact_entries_used * 100. / tot_exact_entries);
  CHECK_EQ(0, tot_histo[0].load()) << "Should have skipped entries with 0 count!";
  size_t tot_histo1 = tot_histo[1].load();
  size_t tot_histo2 = tot_histo[2].load();
  SPLOG("Total kmers with exact counts 1: %ld (%.2f%%): 2: %ld (%.2f%%)", tot_histo1,
        tot_histo1 * 100. / tot_exact_entries, tot_histo2, tot_histo2 * 100. / tot_exact_entries);
  part_used.analyze();
  SPLOG("Per partition min=%.2f%% avg=%.2f%% max=%.2f%% target max=%.2f%%", part_used.min,
        part_used.avg, part_used.max, m_options.max_exact_table_density * 100.);
}

void kmer_counter::extract_exact_counts(
    const std::function<void(extract_iterator, extract_iterator)>& output_f) {
  CHECK_EQ(count_state::EXACT_PASSES_FINISHED, m_count_state);

  SPLOG("Extracting exact counts");
  parallel_for(0, m_options.partitions, [&](size_t part_num) {
    CHECK(m_exact_table[part_num]);
    const auto& et = *m_exact_table[part_num];
    const auto& eto = *m_exact_overflow_table;

    extract_iterator part_begin(et.begin(), &eto);
    extract_iterator part_end(et.end(), &eto);

    output_f(part_begin, part_end);
  });
}

void kmer_counter::close() {
  CHECK_EQ(count_state::EXACT_PASSES_FINISHED, m_count_state);
  m_count_state = count_state::CLOSED;

  SPLOG("kmer_counter: closing");
  m_exact_overflow_table.reset();

  for (unsigned i = 0; i < m_options.partitions; ++i) {
    CHECK(!m_prob_table[i]);
    CHECK(!m_mutable_prob_table[i]);
    m_exact_table[i].reset();
  }
}

kmer_counter::prob_pass_processor::prob_pass_processor(kmer_counter& k)
    : pass_processor(k, 1 /* only one probabilistic pass */) {
  CHECK_EQ(count_state::PROB_PASS, m_kmer_counter.m_count_state);
}

kmer_counter::prob_pass_processor::~prob_pass_processor() {
  CHECK_EQ(count_state::PROB_PASS, m_kmer_counter.m_count_state);
  flush_all();
}

void kmer_counter::prob_pass_processor::flush_part(unsigned part_num, const kmer_t* start,
                                                   const kmer_t* limit) {
  if (start == limit) {
    return;
  }

  auto kmer_size = m_kmer_counter.m_options.kmer_size;
  CHECK(m_kmer_counter.m_mutable_prob_table[part_num]);
  auto& pt = *m_kmer_counter.m_mutable_prob_table[part_num];
  size_t pt_size = pt.size();
  libdivide::divider<uint64_t, libdivide::BRANCHFREE> pt_divider(pt_size);

  kmer_t canon = canonicalize((*start) & k_kmer_mask, kmer_size);
  uint64_t hash = pt_hash_kmer(canon);
  uint64_t pt_pos = hash - pt_divider.perform_divide(hash) * pt_size;
  DCHECK_LT(pt_pos, pt_size);

  kmer_t next_canon = canon;
  uint64_t next_hash = hash;
  uint64_t next_pt_pos = pt_pos;

  do {
    ++start;

    if (start != limit) {
      next_canon = canonicalize((*start) & k_kmer_mask, kmer_size);
      next_hash = pt_hash_kmer(next_canon);
      next_pt_pos = next_hash - pt_divider.perform_divide(next_hash) * pt_size;
      DCHECK_LT(next_pt_pos, pt_size);
      pt.at(next_pt_pos).prefetch_write();
    }

    pt.at(pt_pos).safe_increment();

    canon = next_canon;
    hash = next_hash;
    pt_pos = next_pt_pos;
  } while (start != limit);
}

kmer_counter::exact_pass_processor::exact_pass_processor(kmer_counter& k)
    : pass_processor(k, k.m_exact_passes) {
  CHECK_EQ(count_state::EXACT_PASS, m_kmer_counter.m_count_state);
}

kmer_counter::exact_pass_processor::~exact_pass_processor() {
  CHECK_EQ(count_state::EXACT_PASS, m_kmer_counter.m_count_state);
  flush_all();
}

void kmer_counter::exact_pass_processor::flush_part(unsigned part_num, const kmer_t* start,
                                                    const kmer_t* limit) {
  if (start == limit) {
    return;
  }

  auto kmer_size = m_kmer_counter.m_options.kmer_size;
  CHECK(m_kmer_counter.m_exact_table[part_num]);
  auto& et = *m_kmer_counter.m_exact_table[part_num];
  CHECK(m_kmer_counter.m_exact_overflow_table);
  auto& eto = *m_kmer_counter.m_exact_overflow_table;

  CHECK(m_kmer_counter.m_prob_table[part_num]);
  const auto& pt = *m_kmer_counter.m_prob_table[part_num];
  size_t pt_size = pt.size();
  libdivide::divider<uint64_t, libdivide::BRANCHFREE> pt_divider(pt_size);

  size_t skipped = 0;
  m_kmer_counter.m_tot_exact_kmers.fetch_add(limit - start);

  unsigned min_count = pt.max_value();
  if (m_kmer_counter.m_options.min_count < min_count) {
    min_count = m_kmer_counter.m_options.min_count;
    CHECK_GT(min_count, 0);
  }

  bool flipped;
  kmer_t kmer_and_flags = (*start);
  kmer_t canon = canonicalize(kmer_and_flags & k_kmer_mask, kmer_size, flipped);
  uint64_t hash = pt_hash_kmer(canon);
  uint64_t pt_pos = hash - pt_divider.perform_divide(hash) * pt_size;

  kmer_t next_kmer_and_flags = kmer_and_flags;
  bool next_flipped = flipped;
  kmer_t next_canon = canon;
  uint64_t next_hash = hash;
  uint64_t next_pt_pos = pt_pos;

  do {
    ++start;

    if (start != limit) {
      next_kmer_and_flags = (*start);
      next_canon = canonicalize(next_kmer_and_flags & k_kmer_mask, kmer_size, next_flipped);
      next_hash = pt_hash_kmer(next_canon);
      next_pt_pos = next_hash - pt_divider.perform_divide(next_hash) * pt_size;

      pt.at(next_pt_pos).prefetch_read();
      et.prefetch_write(next_canon);
    }

    if (pt.at(pt_pos) >= min_count) {
      if (et.increment(canon, flipped, kmer_and_flags & k_fwd_flag, kmer_and_flags & k_rev_flag) ==
          exact_count_table_t::max_value()) {
        CHECK_NE(eto.increment(canon, flipped), exact_overflow_count_table_t::max_value())
            << "Overflow on overflow table?";
      }
    } else {
      skipped++;
    }

    kmer_and_flags = next_kmer_and_flags;
    flipped = next_flipped;
    canon = next_canon;
    hash = next_hash;
    pt_pos = next_pt_pos;
  } while (start != limit);
  m_kmer_counter.m_prob_skipped.fetch_add(skipped);
}

void kmer_counter::set_progress_handler(progress_handler_t progress) {
  m_options.progress = progress;
}

}  // namespace build_seqset

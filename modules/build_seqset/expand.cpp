#include "modules/build_seqset/expand.h"

#include <sys/mman.h>

#include "modules/build_seqset/part_counts.h"
#include "modules/io/autostats.h"
#include "modules/io/make_unique.h"
#include "modules/io/parallel.h"
#include "modules/io/stats.h"
#include "modules/io/track_mem.h"

namespace build_seqset {

namespace {

time_t g_last_log_state = 0;

template <typename it_type>
void skip_dups(it_type& cur, it_type& end, const dna_slice& min_next_entry) {
  if (cur == end) {
    return;
  }
  it_type next = cur;
  ++next;
  while (next != end) {
    auto cmp = cur->compare_to(*next);
    switch (cmp) {
      case dna_compare_result::FIRST_IS_PREFIX:
      case dna_compare_result::EQUAL:
        ++cur;
        ++next;
        continue;
      default:
        return;
    }
  }
  if (min_next_entry.size() > 0) {
    auto cmp = cur->compare_to(min_next_entry);
    switch (cmp) {
      case dna_compare_result::FIRST_IS_PREFIX:
      case dna_compare_result::EQUAL:
        ++cur;
        break;
      default:
        return;
    }
  }
}

struct part_expander_stats : public autostats_base {
  DECLARE_AUTOSTATS(part_expander_stats,          //
                    ((COUNTER, new_entries))      //
                    ((COUNTER, sorted_entries))   //
                    ((COUNTER, sorted_output))    //
                    ((COUNTER, output_entries))   //
                    ((COUNTER, expanded_output))  //
                    ((COUNTER, prefetch_msecs))   //
                    ((COUNTER, sort_msecs))       //
                    ((COUNTER, dedup_msecs))      //
                    );

  std::string as_string() {
    std::stringstream result;
    write_to_stream(result);
    return result.str();
  }
};

size_t to_msecs(std::chrono::duration<double> d) {
  return size_t(d.count() * 1000);
}

class part_expander {
 public:
  part_expander(part_repo& entries, part_counts* counts, part_repo::partition_ref& sorted_part,
                part_repo::partition_ref& new_part, kmer_t part_id, unsigned expand_count,
                unsigned expand_stride, const std::string& sorted_out_pass,
                const std::string& expanded_out_pass)
      : m_entries(entries),
        m_counts(counts),
        m_sorted_part(sorted_part),
        m_new_part(new_part),
        m_part_id(part_id),
        m_expand_count(expand_count),
        m_expand_stride(expand_stride),
        m_sorted_out_pass(sorted_out_pass),
        m_expanded_out_pass(expanded_out_pass) {
    CHECK_EQ(m_new_part.prefix, m_sorted_part.prefix);
    CHECK(m_new_part.main->repo().begin() == m_sorted_part.main->repo().begin());

    m_empty = m_new_part.main->begin() == m_new_part.main->end() &&
              m_sorted_part.main->begin() == m_sorted_part.main->end();

    if (m_empty) {
      m_new_part.reset();
      m_sorted_part.reset();
      return;
    }

    m_num_new_entries = m_new_part.main->end() - m_new_part.main->begin();
    m_stats.new_entries += m_num_new_entries;
    m_num_sorted_entries = m_sorted_part.main->end() - m_sorted_part.main->begin();
    m_stats.sorted_entries += m_num_sorted_entries;
  }

  bool empty() const {
    return m_empty;
  }

  size_t new_entries_memory() const {
    return sizeof(seq_repository::entry_data) * m_num_new_entries;
  }
  size_t sorted_entries_memory() const {
    return sizeof(seq_repository::entry_data) * m_num_sorted_entries;
  }

  size_t sort_memory_needed() const { return new_entries_memory(); }

  size_t max_memory_needed() const { return std::max(sort_memory_needed(), dedup_memory_needed()); }

  size_t dedup_memory_needed() const { return new_entries_memory() + sorted_entries_memory(); }

  void do_output(parallel_state& st) {
    set_state("output");
    st.unreserve_memory(max_memory_needed() - dedup_memory_needed());
    dedup_and_output();
    CHECK_EQ(m_new_data.size(), 0) << "Memory should be free after end of expand";
    st.unreserve_memory(dedup_memory_needed());
    set_state("done");
  }

  void set_state(const std::string& desc) {
    constexpr int k_set_state_seconds = 60;
    static std::mutex set_state_mu;
    std::lock_guard<std::mutex> l(set_state_mu);
    std::string part_name = m_new_part.prefix.as_string();

    time_t now = time(0);
    auto& st = g_part_states[part_name];
    if (st.second.empty()) {
      st.first = now;
    }
    if (st.second == desc) {
      CHECK_NE(desc, "done");
      return;
    }
    st.second = desc;

    if (g_last_log_state + k_set_state_seconds <= now) {
      if (g_last_log_state) {
        log_state();
      }
      g_last_log_state = now;
    }
    if (desc == "done") {
      g_part_states.erase(part_name);
      if (g_part_states.empty()) {
        g_last_log_state = 0;
      }
    } else {
      st.first = now;
    }
  }

  void log_state() {
    time_t now = time(0);

    std::stringstream out;
    for (const auto& ps : g_part_states) {
      const auto& part_name = ps.first;
      time_t start_time = ps.second.first;
      const std::string& state = ps.second.second;
      int duration = now - start_time;

      out << " " << part_name << ":" << state << "(" << duration << "s)";
    }
    SPLOG("%lu parts in progress:%s", g_part_states.size(), out.str().c_str());
  }

  void do_sort() {
    CHECK(!m_empty);
    CHECK(m_new_data_begin);

    set_state("sort");

    auto sort_start_time = std::chrono::steady_clock::now();

    parallel_for(0, m_section_starts.size(), [&](size_t i) {
      m_sorted_part.main->sort_entry_data(m_section_starts[i], m_section_limits[i]);
    });

    auto sort_end_time = std::chrono::steady_clock::now();

    m_stats.sort_msecs += to_msecs(sort_end_time - sort_start_time);
  }

  part_expander_stats stats() const { return m_stats; }

  void do_prefetch() {
    if (m_empty) {
      return;
    }

    auto prefetch_start_time = std::chrono::steady_clock::now();

    set_state("prefetch");

    m_new_data = mutable_membuf(new owned_membuf(
        m_num_new_entries * sizeof(seq_repository::entry_data), "build_seqset_expand_sort"));
    m_new_data_begin = (seq_repository::entry_data*)m_new_data.mutable_data();
    m_new_data_end = m_new_data_begin + m_num_new_entries;

    CHECK_EQ(reinterpret_cast<const char*>(m_new_data_end),
             m_new_data.mutable_data() + m_new_data.size());

    if (m_counts) {
      auto count_index_bounds = m_counts->seq_to_index_range(m_new_part.prefix);

      CHECK_NE(count_index_bounds.first, count_index_bounds.second);

      m_section_starts.reserve(count_index_bounds.second - count_index_bounds.first);

      size_t cur_idx = 0;
      for (size_t i = count_index_bounds.first; i != count_index_bounds.second; ++i) {
        m_section_starts.push_back(m_new_data_begin + cur_idx);
        cur_idx += m_counts->counts()[i];
      }
      if (cur_idx != size_t(m_new_data_end - m_new_data_begin)) {
        auto it = m_new_part.main->data_begin();
        auto prefetch_end = m_new_part.main->data_end();
        while (it != prefetch_end) {
          ++it;
        }
      }
      CHECK_EQ(cur_idx, (m_new_data_end - m_new_data_begin));
      CHECK_GE(m_counts->bases(), m_entries.partition_depth());
      m_section_bases = m_counts->bases() - m_entries.partition_depth();
    } else {
      m_section_starts.push_back(m_new_data_begin);
      m_section_bases = 0;
    }

    CHECK_EQ(m_section_starts.size(), 1UL << (m_section_bases * 2));

    auto it = m_new_part.main->data_begin();
    auto prefetch_end = m_new_part.main->data_end();

    madvise(const_cast<void*>(reinterpret_cast<const void*>(it)),
            reinterpret_cast<const char*>(prefetch_end) - reinterpret_cast<const char*>(it),
            MADV_SEQUENTIAL);

    m_section_limits = m_section_starts;

    while (it != prefetch_end) {
      const auto& e = *it;
      auto dna_it =
          dna_const_iterator(e.raw_inline_bases(), m_entries.partition_depth(), 0 /* rev comp */);

      size_t section = 0;
      for (size_t b_idx = 0; b_idx != m_section_bases; ++b_idx) {
        section <<= 2;
        section |= int(dna_base(*dna_it));
        ++dna_it;
      }

      *m_section_limits[section]++ = e;
      ++it;
    }

    CHECK_EQ(m_section_starts[0], m_new_data_begin);
    for (size_t i = 0; i != m_section_starts.size(); ++i) {
      if (i + 1 == m_section_starts.size()) {
        CHECK_EQ(m_section_limits[i], m_new_data_end);
      } else {
        CHECK_EQ(m_section_limits[i], m_section_starts[i + 1]);
      }
    }

    m_new_part.reset();
    auto prefetch_end_time = std::chrono::steady_clock::now();
    m_stats.prefetch_msecs += to_msecs(prefetch_end_time - prefetch_start_time);
  }
 public:

  void dedup_and_output() {
    auto dedup_start_time = std::chrono::steady_clock::now();

    auto new_begin = seq_repository::iterator(m_new_data_begin, m_sorted_part.main->repo());
    auto new_end = new_begin + m_num_new_entries;
    auto new_cur = new_begin;

    auto sorted_begin = m_sorted_part.main->begin();
    auto sorted_end = m_sorted_part.main->end();
    auto sorted_cur = sorted_begin;

    madvise(const_cast<void*>(reinterpret_cast<const void*>(m_sorted_part.main->data_begin())),
            reinterpret_cast<const char*>(m_sorted_part.main->data_end()) -
                reinterpret_cast<const char*>(m_sorted_part.main->data_begin()),
            MADV_SEQUENTIAL);

    auto output_builder = m_entries.open_ref_builder(m_part_id, m_sorted_out_pass);
    std::function<void(const seq_repository::entry_base& e)> output =
        [&](const seq_repository::entry_base& e) {
          m_stats.output_entries++;
          output_builder->write_entry_unlocked(e);
        };

    std::function<void(const seq_repository::entry_base& e)> output_and_expand =
        [&](const seq_repository::entry_base& e) {
          m_stats.output_entries++;
          output_builder->write_entry_unlocked(e);
          if (m_expand_count) {
            m_stats.expanded_output +=
                m_entries.write_with_expansions(e.pop_front(), m_expand_stride, m_expand_count);
          }
        };

    dna_sequence min_next_entry;
    if (m_sorted_part.next_entry.size() > 0) {
      if (m_new_part.next_entry.size() > 0 && m_new_part.next_entry < m_sorted_part.next_entry) {
        min_next_entry = m_new_part.next_entry;
      } else {
        min_next_entry = m_sorted_part.next_entry;
      }
    } else {
      min_next_entry = m_new_part.next_entry;
    }

    skip_dups(new_cur, new_end, min_next_entry);
    while (new_cur != new_end || sorted_cur != sorted_end) {
      dna_compare_result cmp;

      if (new_cur == new_end) {
        if (min_next_entry.size() > 0) {
          cmp = sorted_cur->compare_to(min_next_entry);
          CHECK_NE(cmp, dna_compare_result::SECOND_IS_LESS);
          CHECK_NE(cmp, dna_compare_result::SECOND_IS_PREFIX);
          CHECK_NE(cmp, dna_compare_result::EQUAL);
        } else {
          cmp = dna_compare_result::FIRST_IS_LESS;
        }
      } else if (sorted_cur == sorted_end) {
        cmp = dna_compare_result::SECOND_IS_LESS;
      } else {
        cmp = sorted_cur->compare_to(*new_cur);
      }
      // First is sorted, Second is new.
      switch (cmp) {
        case dna_compare_result::FIRST_IS_PREFIX:
          CHECK(sorted_cur != sorted_end);
          ++sorted_cur;
          break;
        case dna_compare_result::FIRST_IS_LESS:
          CHECK(sorted_cur != sorted_end);
          output(*sorted_cur);
          ++sorted_cur;
          break;
        case dna_compare_result::SECOND_IS_LESS:
          CHECK(new_cur != new_end);
          output_and_expand(*new_cur);
          ++new_cur;
          skip_dups(new_cur, new_end, min_next_entry);
          break;
        case dna_compare_result::SECOND_IS_PREFIX:
        case dna_compare_result::EQUAL:
          // Duplicate; ignore
          CHECK(new_cur != new_end);
          ++new_cur;
          skip_dups(new_cur, new_end, min_next_entry);
          break;
      }
    }

    m_new_data_begin = nullptr;
    m_new_data_end = nullptr;
    m_new_data = mutable_membuf();
    m_sorted_part.reset();
    auto dedup_end_time = std::chrono::steady_clock::now();
    m_stats.dedup_msecs += to_msecs(dedup_end_time - dedup_start_time);
  }

  part_repo& m_entries;
  part_counts* m_counts;
  part_repo::partition_ref& m_sorted_part;
  part_repo::partition_ref& m_new_part;

  kmer_t m_part_id;
  unsigned m_expand_count;
  unsigned m_expand_stride;
  std::string m_sorted_out_pass;
  std::string m_expanded_out_pass;

  part_expander_stats m_stats;
  size_t m_num_new_entries = 0;
  size_t m_num_sorted_entries = 0;

  bool m_empty = false;

  size_t m_section_bases = 0;
  std::vector<seq_repository::entry_data*> m_section_starts;
  std::vector<seq_repository::entry_data*> m_section_limits;

  // New unsorted, from prefetch:
  mutable_membuf m_new_data;
  seq_repository::entry_data* m_new_data_begin = nullptr;
  seq_repository::entry_data* m_new_data_end = nullptr;

  // Current state of each partition
  using part_states_t = std::map<std::string /* partition */,
                                 std::pair<time_t /* start time */, std::string /* state */>>;
  static part_states_t g_part_states;
};

part_expander::part_states_t part_expander::g_part_states;

}  // namespace

size_t expander::sort_and_dedup(const std::string& already_sorted_pass,
                                const std::string& new_entries_pass,
                                const std::string& sorted_out_pass,
                                const std::string& expanded_out_pass, unsigned expand_stride,
                                unsigned expand_count, progress_handler_t progress) {
  progress(0);

  size_t seq_repo_bases = m_entries.repo_slice().size();
  size_t seq_repo_size = seq_repo_bases / 4;
  size_t memory_bytes = get_maximum_mem_bytes();
  if (seq_repo_size > memory_bytes) {
    SPLOG(
        "WARNING: Sequence repo (%.2f MB) takes up more than RAM available "
        "(%.2f MB).  Pretending it's smaller.",
        seq_repo_size / 1024. / 1024, memory_bytes / 1024. / 1024.);
    seq_repo_size = memory_bytes;
  }

  memory_bytes -= seq_repo_size;

  SPLOG(
      "Sorting and deduping \"%s\" + \"%s\" -> \"%s\" + \"%s\", reserving "
      "%.2f MB RAM for sequence repo, limiting sort/dedup to %.2f MB RAM",
      already_sorted_pass.c_str(), new_entries_pass.c_str(), sorted_out_pass.c_str(),
      expanded_out_pass.c_str(), seq_repo_size / 1024. / 1024, memory_bytes / 1024. / 1024);
  track_mem::reset_stats();

  m_part_counts = m_entries.release_part_counts(new_entries_pass);
  if (!m_part_counts) {
    SPLOG("WARNING: Part counts not available; this stage may be slow");
  }

  auto sorted_parts = m_entries.partitions(already_sorted_pass, false /* no expansions nedeed */,
                                           !m_keep_tmp /* delete on close */);
  auto new_parts = m_entries.partitions(new_entries_pass, false /* no expansions needed */,
                                        !m_keep_tmp /* delete on close */);
  if (expand_count) {
    CHECK_NE(expanded_out_pass, "");
    m_entries.open_write_pass(expanded_out_pass);
  }

  size_t num_parts = new_parts.size();
  std::vector<std::unique_ptr<part_expander>> expanders(num_parts);

  for (kmer_t part_id = 0; part_id < num_parts; ++part_id) {
    expanders[part_id] = make_unique<part_expander>(
        m_entries, m_part_counts.get(), sorted_parts[part_id], new_parts[part_id], part_id,
        expand_count, expand_stride, sorted_out_pass, expanded_out_pass);
  }

  if (true) {
    auto less_memory_first = [](const std::unique_ptr<part_expander>& lhs,
                                const std::unique_ptr<part_expander>& rhs) {
      return lhs->max_memory_needed() < rhs->max_memory_needed();
    };
    auto more_memory_first = [](const std::unique_ptr<part_expander>& lhs,
                                const std::unique_ptr<part_expander>& rhs) {
      return lhs->max_memory_needed() > rhs->max_memory_needed();
    };

    unsigned first_small = get_thread_count() / 2;
    if (first_small < 4) {
      first_small = 4;
    }

    if (expanders.size() > first_small) {
      auto it = expanders.begin() + first_small;
      std::nth_element(expanders.begin(), it, expanders.end());
      std::sort(expanders.begin(), it, less_memory_first);
      std::sort(it, expanders.end(), more_memory_first);
    } else {
      std::sort(expanders.begin(), expanders.end(), less_memory_first);
    }
  } else {
    // Sort the smallest ones first so we can ramp up quickly.
    std::sort(
        expanders.begin(), expanders.end(),
        [](const std::unique_ptr<part_expander>& lhs, const std::unique_ptr<part_expander>& rhs) {
          return lhs->max_memory_needed() < rhs->max_memory_needed();
        });

    CHECK(!expanders.empty());
    // Distribute the bigger ones throughout the workload with the biggest ones towards the
    // beginning.
    unsigned biggest = expanders.size() - 1;
    for (unsigned i = expanders.size() / 2; i < biggest; i += get_thread_count(), --biggest) {
      std::swap(expanders[i], expanders[biggest]);
    }
  }

  CHECK_EQ(num_parts, sorted_parts.size());
  CHECK(!expanders.empty());

  std::mutex mu;
  std::deque<thread_pool::work_t> deferred_worklist;
  int cur_prio = 0;
  bool more_prefetch_queued = false;
  auto queue_more_if_needed = [&](bool if_not_already_queued = false) {
    std::lock_guard<std::mutex> l(mu);
    if (deferred_worklist.empty()) {
      return;
    }
    if (more_prefetch_queued) {
      // Only have a single prefetch queued to be started at once.
      return;
    }
    more_prefetch_queued = true;
    --cur_prio;
    thread_pool::work_t w = std::move(deferred_worklist.front());
    w.progress_part = 1. - 1. / deferred_worklist.size();
    deferred_worklist.pop_front();
    parallel_pool().add_work_async(std::move(w), cur_prio);
  };
  for (auto& e : expanders) {
    if (e->empty()) {
      continue;
    }
    thread_pool::work_t work{
        [&e, queue_more_if_needed, &mu, &more_prefetch_queued](parallel_state& st) {
          {
            std::lock_guard<std::mutex> l(mu);
            // We're executing a queued prefetch.
            more_prefetch_queued = false;
          }
          e->do_prefetch();
          queue_more_if_needed();
          e->do_sort();
          e->do_output(st);
          queue_more_if_needed();
        }};
    work.reserve_memory = e->max_memory_needed();
    if (work.reserve_memory > memory_bytes) {
      SPLOG("WARNING: Increasing max memory from %ld to %ld to accomodate large part", memory_bytes,
            work.reserve_memory);
      memory_bytes = work.reserve_memory;
    }
    deferred_worklist.emplace_back(std::move(work));
  }

  parallel_pool().set_memory_limit(memory_bytes);
  parallel_pool().execute_worklist(
      std::vector<thread_pool::work_t>{
          {[queue_more_if_needed](parallel_state&) { queue_more_if_needed(); }}},
      progress);
  CHECK(deferred_worklist.empty());

  part_expander_stats tot_stats;
  for (unsigned i = 0; i < num_parts; ++i) {
    tot_stats += expanders[i]->stats();
  }
  expanders.clear();
  m_entries.flush();
  // Don't reuse part counts for the next pass.
  m_part_counts.reset();

  SPLOG("Expand stats: %s", tot_stats.as_string().c_str());

  CHECK_GE(tot_stats.new_entries + tot_stats.sorted_entries, tot_stats.output_entries);
  size_t tot_dedupped = tot_stats.new_entries + tot_stats.sorted_entries - tot_stats.output_entries;

  return tot_dedupped;
}

size_t expander::expand(const std::string& input_pass,
                        const std::string& expanded_pass, unsigned stride,
                        unsigned count, progress_handler_t progress) {
  std::atomic<size_t> tot_expanded{0};
  std::atomic<size_t> tot_entries{0};

  SPLOG("Expanding with stride=%d, count=%d", stride, count);
  track_mem::reset_stats();

  m_entries.open_write_pass(expanded_pass);
  m_entries.for_each_partition(  //
      input_pass,
      [&](const part_repo::partition_ref& part) {
        size_t chunk_expanded = 0;

        dna_base_array<seq_repository::popped_iterator> bcur;
        dna_base_array<seq_repository::popped_iterator> bend;

        for (dna_base b : dna_bases()) {
          bcur[b] = part.pushed[b].first.pop_front();
          bend[b] = part.pushed[b].second.pop_front();
        }
        auto it_end = part.main->end();
        auto it = part.main->begin();
        tot_entries.fetch_add(it_end - it);
        while (it != it_end) {
          for (dna_base b : dna_bases()) {
            seq_repository::popped_iterator& base_cur = bcur[b];
            const seq_repository::popped_iterator& base_end = bend[b];

            bool done_advancing = false;
            while (!done_advancing) {
              if (base_cur == base_end) {
                done_advancing = true;
                break;
              }

              auto cmp = base_cur->compare_to(*it);
              switch (cmp) {
                case dna_compare_result::FIRST_IS_LESS:
                  chunk_expanded +=
                      m_entries.write_with_expansions(*base_cur, stride, count);
                  done_advancing = false;
                  break;
                case dna_compare_result::FIRST_IS_PREFIX:
                case dna_compare_result::EQUAL: {
                  ++base_cur;
                  done_advancing = true;
                  break;
                }
                case dna_compare_result::SECOND_IS_LESS:
                case dna_compare_result::SECOND_IS_PREFIX:
                  done_advancing = true;
                  break;
              }
              if (!done_advancing) {
                ++base_cur;
              }
            }
          }
          ++it;
        }
        for (dna_base b : dna_bases()) {
          seq_repository::popped_iterator& base_cur = bcur[b];
          const seq_repository::popped_iterator& base_end = bend[b];

          while (base_cur != base_end) {
            chunk_expanded +=
                m_entries.write_with_expansions(*base_cur, stride, count);
            ++base_cur;
          }
        }
        tot_expanded += chunk_expanded;
      },
      progress);

  m_entries.flush();
  SPLOG("Expand pass completed with %ld (%.2f%%) new entries",
        tot_expanded.load(), tot_expanded.load() * 100. / (tot_entries.load()));
  return tot_expanded.load();
}

}  // namespace build_seqset

#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/function_output_iterator.hpp>

#include "modules/io/file_io.h"
#include "modules/variants/read_cov.h"

#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"

#include "absl/container/btree_map.h"

namespace variants {

constexpr bool k_cov_debug = false;
constexpr bool k_show_stats = false;

struct read_cov::result {
  assembly_ptr a;
  size_t idx;
  aoffset_t seq_size;
  read_coverage_set reads;
  read_cov* cov = nullptr;
};

struct read_cov::keep_result {
  aoffset_t keep_until = std::numeric_limits<aoffset_t>::min();
  result* r = nullptr;

  bool operator<(const keep_result& rhs) const { return r < rhs.r; }
};

struct read_cov::result_offset {
  aoffset_t end_offset;
  result* r = nullptr;

  aoffset_t start_offset() const {
    CHECK(r);
    return end_offset - aoffset_t(r->seq_size);
  }
  bool operator==(const result_offset& rhs) const {
    return end_offset == rhs.end_offset && r == rhs.r;
  }
  bool operator<(const result_offset& rhs) const {
    CHECK(r);
    CHECK(rhs.r);
    if (end_offset != rhs.end_offset) {
      return end_offset < rhs.end_offset;
    }
    if (r != rhs.r) {
      return r->idx < rhs.r->idx;
    }
    return false;
  }
};

class read_cov::cov_pg {
 public:
  using keep_results_t = std::vector<keep_result>;
  cov_pg(read_cov* cov, const seqset_range& r);
  ~cov_pg() { sub_stats(); }

  std::unique_ptr<cov_pg> split();
  void join(std::unique_ptr<cov_pg> other);

  void add_result(std::unique_ptr<result> r);
  roff_index_t add_result_offset(result_offset roff);

  void set_max_size(size_t new_size) { m_max_size = new_size; }
  void adjust_pos(aoffset_t adjust);

  void add_stats();
  void sub_stats();

  void compact();
  void compact_internal();
  const keep_results_t& keep_results() const { return m_keep_results; }

  void save_path_priorities();

 private:
  class path_info {
   public:
    static constexpr size_t k_prio_not_used = std::numeric_limits<size_t>::max();
    static constexpr size_t k_init_prio = 1;

    path_info() = default;
    path_info(const path_info&) = delete;
    path_info(roff_index_t* begin, roff_index_t* end, size_t prio)
        : m_begin(begin), m_end(end), m_path_priority(prio) {
      CHECK(m_begin);
      CHECK(m_end);
      CHECK_GE(m_end, m_begin);
    }
    path_info(path_info&& rhs)
        : m_begin(rhs.m_begin), m_end(rhs.m_end), m_path_priority(rhs.m_path_priority) {
      rhs.m_begin = nullptr;
      rhs.m_end = nullptr;
      rhs.m_path_priority = k_prio_not_used;
      CHECK_NE(m_path_priority, k_prio_not_used);
    }
    path_info& operator=(const path_info&) = delete;
    path_info& operator=(path_info&& rhs) {
      if (this != &rhs) {
        m_begin = rhs.m_begin;
        m_end = rhs.m_end;
        m_path_priority = rhs.m_path_priority;
        rhs.m_begin = nullptr;
        rhs.m_end = nullptr;
        rhs.m_path_priority = k_prio_not_used;
        CHECK_NE(m_path_priority, k_prio_not_used);
      }
      return *this;
    }

    size_t size() const {
      DCHECK(m_begin);
      DCHECK(m_end);
      DCHECK_GE(m_end, m_begin);
      return m_end - m_begin;
    }
    bool empty() const { return m_begin == m_end; }

    roff_index_t* begin() { return m_begin; }
    roff_index_t* end() { return m_end; }
    const roff_index_t* begin() const { return m_begin; }
    const roff_index_t* end() const { return m_end; }

    void set_end(roff_index_t* new_end) {
      CHECK_GE(new_end, m_begin);
      CHECK_LE(new_end, m_end);
      m_end = new_end;
    }

    void increase_path_priority(size_t new_priority) {
      CHECK_NE(m_path_priority, k_prio_not_used);
      m_path_priority = std::max(m_path_priority, new_priority);
    }
    size_t path_priority() const {
      CHECK_NE(m_path_priority, k_prio_not_used);
      return m_path_priority;
    }

   private:
    roff_index_t* m_begin = nullptr;
    roff_index_t* m_end = nullptr;

    // This path would have been trimmed if assemble_options::max_coverage_paths was less than this
    // value:
    size_t m_path_priority = k_prio_not_used;
  };

  struct pending_reads_key {
    size_t roff_idx;
    aoffset_t offset;
    aoffset_t read_len;

    template <typename H>
    friend H AbslHashValue(H h, const pending_reads_key& k) {
      return H::combine(std::move(h), k.offset, k.read_len);
    }

    bool operator==(const pending_reads_key& rhs) const {
      return roff_idx == rhs.roff_idx && offset == rhs.offset && read_len == rhs.read_len;
    }
  };

  struct pending_reads_val {
    read_id_set read_ids;
    size_t max_paths = 0;
  };

  using paths_elem_t = std::pair<seqset_range, path_info>;
  using paths_t = std::vector<paths_elem_t>;
  using roff_translate_table_t = std::vector<roff_index_t>;

  // Returns the total number of roffrefs, including unmerged.
  size_t count_roffrefs() const;
  // num_roffrefs is the result from count_roffrefs().
  void copy_paths_from(size_t num_roffrefs, const paths_t& old_paths,
                       const std::vector<roff_index_t>& unmerged_roffrefs);
  void add_base(dna_base b);
  void add_path(const seqset_range& r, path_info pi,
                const roff_translate_table_t* translate = nullptr,
                const std::vector<roff_index_t>* unmerged_roffrefs = nullptr)
      __attribute__((noinline));
  void remove_expired_roffs(path_info& pi, aoffset_t r_size);
  void flush_pending_reads();
  void flush_old_results();
  void trim_paths_to_max();
  void make_roff_translates(roffs_t roffs1, roffs_t roffs2, roff_translate_table_t& rt1,
                            roff_translate_table_t& rt2);
  static keep_results_t join_keep_results(keep_results_t, keep_results_t);

  path_info alloc_pi(size_t nelem, size_t prio);
  path_info alloc_big_pi(size_t nelem, size_t prio);
  void unalloc_pi(roff_index_t* begin, roff_index_t* end);

  cov_pg() = default;

  read_cov* m_cov = nullptr;
  aoffset_t m_pos = 0;
  roffs_t m_result_offsets;
  std::vector<roff_index_t> m_unmerged_roffrefs;
  keep_results_t m_keep_results;
  bool m_need_compact = false;

  // Sorted by seqset range.
  struct paths_comparer;
  paths_t m_paths;

  size_t m_max_size = 0;

  // Allocations for roff indexes that are of sizes other than
  // k_pi_alloc_chunk.
  std::vector<std::unique_ptr<roff_index_t[]>> m_big_alloc_pi;
  // Sum of the sizes of the blocks in m_big_alloc_pi.
  size_t m_big_alloc_pi_size = 0;

  // Used blocks of roff indexes, in chunks of k_pi_alloc_chunk elements.
  std::vector<std::unique_ptr<roff_index_t[]>> m_alloc_pi;
  uint32_t* m_alloc_pi_start = nullptr;
  size_t m_alloc_pi_free = 0;

  absl::flat_hash_map<pending_reads_key, pending_reads_val> m_pending_reads;
};

constexpr size_t read_cov::k_pi_alloc_chunk;
constexpr size_t read_cov::cov_pg::path_info::k_prio_not_used;
constexpr size_t read_cov::cov_pg::path_info::k_init_prio;

void read_cov::cov_pg::add_stats() {
  if (!k_show_stats) {
    return;
  }
  ++m_cov->m_tot_pg;
  m_cov->m_tot_results += m_keep_results.size();
  m_cov->m_tot_paths += m_paths.size();
  m_cov->m_tot_roffs += m_result_offsets.size();
  m_cov->m_tot_roff_refs += count_roffrefs();
  m_cov->m_tot_pialloc_blocks += m_alloc_pi.size();
  m_cov->m_tot_big_pialloc_blocks += m_big_alloc_pi.size();
  m_cov->m_tot_big_pialloc_size += m_big_alloc_pi_size;
}

void read_cov::cov_pg::sub_stats() {
  if (!k_show_stats) {
    return;
  }
  CHECK_GE(m_cov->m_tot_pg, 1);
  --m_cov->m_tot_pg;

  size_t tot_results = m_keep_results.size();
  CHECK_GE(m_cov->m_tot_results, tot_results);
  m_cov->m_tot_results -= tot_results;

  size_t tot_paths = m_paths.size();
  CHECK_GE(m_cov->m_tot_paths, tot_paths);
  m_cov->m_tot_paths -= tot_paths;

  size_t tot_roffs = m_result_offsets.size();
  CHECK_GE(m_cov->m_tot_roffs, tot_roffs);
  m_cov->m_tot_roffs -= tot_roffs;

  size_t tot_roff_refs = count_roffrefs();
  CHECK_GE(m_cov->m_tot_roff_refs, tot_roff_refs);
  m_cov->m_tot_roff_refs -= tot_roff_refs;

  CHECK_GE(m_cov->m_tot_pialloc_blocks, m_alloc_pi.size());
  m_cov->m_tot_pialloc_blocks -= m_alloc_pi.size();
  CHECK_GE(m_cov->m_tot_big_pialloc_blocks, m_big_alloc_pi.size());
  m_cov->m_tot_big_pialloc_blocks -= m_big_alloc_pi.size();
  CHECK_GE(m_cov->m_tot_big_pialloc_size, m_big_alloc_pi_size);
  m_cov->m_tot_big_pialloc_size -= m_big_alloc_pi_size;

  m_cov->on_add_stats();
}

void read_cov::on_add_stats() {
  time_t now = time(0);
  if (now > (m_last_stats_report + 60)) {
    display_stats_report();
  }
}

void read_cov::display_stats_report() {
  size_t keep_results_mb = m_tot_results * sizeof(keep_result) / (1024 * 1024);
  size_t roffs_mb = m_tot_roffs * sizeof(result_offset) / (1024 * 1024);
  size_t roff_refs_mb = m_tot_roff_refs * sizeof(roff_index_t) / (1024 * 1024);
  size_t tot_mb = keep_results_mb + roffs_mb + roff_refs_mb;
  size_t free_pi_mb =
      m_free_pi_allocs.size() * k_pi_alloc_chunk * sizeof(roff_index_t) / (1024 * 1024);
  size_t big_pi_mb = m_tot_big_pialloc_size * sizeof(roff_index_t) / (1024 * 1024);

  std::cerr << "read_cov stats " << tot_mb << " MB cur=" << m_cur_offset
            << " assemblies seen=" << m_assemblies_seen << " pgs=" << m_tot_pg
            << " results=" << m_tot_results << "  (" << keep_results_mb << " MB) roffs "
            << m_tot_roffs << " (" << roffs_mb << " MB) paths=" << m_tot_paths
            << " roff refs=" << m_tot_roff_refs << " (" << roff_refs_mb
            << " MB), pi blocks= " << m_tot_pialloc_blocks << " big=" << m_tot_big_pialloc_blocks
            << " (" << big_pi_mb << " MB), free= " << m_free_pi_allocs.size() << " (" << free_pi_mb
            << " MB)\n";
  m_last_stats_report = time(0);
}

struct read_cov::cov_pg::paths_comparer {
  bool operator()(const std::pair<seqset_range, path_info>& lhs,
                  const std::pair<seqset_range, path_info>& rhs) const {
    if (lhs.first.begin() != rhs.first.begin()) {
      return lhs.first.begin() < rhs.first.begin();
    }
    if (lhs.first.size() != rhs.first.size()) {
      return lhs.first.size() < rhs.first.size();
    }
    return false;
  }
};

void read_cov::translate_uint32s(uint32_t* start, uint32_t* end,
                                 const std::vector<uint32_t>& translate_table) {
  if (start == end) {
    return;
  }
  CHECK(start);
  CHECK_LE(start, end);
  constexpr bool k_fast_translate = true;
  // The following two methods of implementation are equivalent in
  // functionality, but loading and storing 64 bits at a time
  // benchmarks at about twice as fast.  Thanks to
  // https://travisdowns.github.io/blog/2019/06/11/speed-limits.html
  // for ideas here.
  if (k_fast_translate) {
    size_t nelem64 = (end - start) / 2;
    uint32_t* itr64 = start;
    uint32_t* itr64_end = itr64 + nelem64 * 2;

    while (itr64 != itr64_end) {
      uint64_t idx64;
      memcpy(&idx64, itr64, sizeof(idx64));
      uint32_t idx1 = idx64 >> 32;
      uint32_t idx2 = idx64;
      DCHECK_LT(idx1, translate_table.size());
      uint32_t translated1 = translate_table[idx1];
      DCHECK_LT(idx2, translate_table.size());
      uint32_t translated2 = translate_table[idx2];
      uint64_t together = (uint64_t(translated1) << 32) | uint64_t(translated2);
      memcpy(itr64, &together, sizeof(together));

      itr64 += 2;
    }

    start = itr64_end;
  }
  while (start != end) {
    uint32_t idx = *start;
    DCHECK_LT(idx, translate_table.size());
    uint32_t translated = translate_table[idx];
    *start = translated;
    ++start;
  }
}

void read_cov::cov_pg::add_path(const seqset_range& r, path_info new_pi,
                                const roff_translate_table_t* translate,
                                const std::vector<roff_index_t>* unmerged) {
  if (translate) {
    translate_uint32s(new_pi.begin(), new_pi.end(), *translate);
  }

  if (unmerged && !unmerged->empty()) {
    path_info dest_pi = alloc_pi(new_pi.size() + unmerged->size(), 1);
    m_need_compact = true;
    auto it = std::copy(new_pi.begin(), new_pi.end(), dest_pi.begin());
    it = std::copy(unmerged->begin(), unmerged->end(), it);
    CHECK_EQ(it, dest_pi.end());
    new_pi = std::move(dest_pi);
  }

  if (m_paths.empty() || m_paths.back().first != r) {
    m_paths.emplace_back(r, std::move(new_pi));
    return;
  }

  auto& old_pi = m_paths.back().second;
  CHECK(m_paths.back().first == r) << r.sequence() << " vs " << m_paths.back().first.sequence();
  if (old_pi.empty()) {
    old_pi = std::move(new_pi);
  } else if (!new_pi.empty()) {
    path_info dest_pi = alloc_pi(old_pi.size() + new_pi.size(),
                                 std::min(old_pi.path_priority(), new_pi.path_priority()));

    m_need_compact = true;
    auto last_dest =
        std::set_union(old_pi.begin(), old_pi.end(), new_pi.begin(), new_pi.end(), dest_pi.begin());
    unalloc_pi(last_dest, dest_pi.end());
    dest_pi.set_end(last_dest);
    old_pi = std::move(dest_pi);
    DCHECK(std::is_sorted(old_pi.begin(), old_pi.end()));
  }
}

read_cov::~read_cov() {
  flush();
  CHECK_EQ(0, m_tot_pg);
  CHECK_EQ(0, m_tot_results);
  CHECK_EQ(0, m_tot_roffs);
  CHECK_EQ(0, m_tot_paths);
  CHECK_EQ(0, m_tot_roff_refs);
  CHECK_EQ(0, m_tot_pialloc_blocks);
  CHECK_EQ(0, m_tot_big_pialloc_blocks);
  CHECK_EQ(0, m_tot_big_pialloc_size);
}

read_cov::cov_pg::cov_pg(read_cov* cov, const seqset_range& r) : m_cov(cov) {
  m_paths.emplace_back(r, alloc_pi(0, path_info::k_init_prio));
  add_stats();
}

read_cov::roff_index_t read_cov::cov_pg::add_result_offset(result_offset roff) {
  if (!m_result_offsets.empty()) {
    DCHECK_LT(m_result_offsets.back(), roff);
  }
  m_result_offsets.push_back(std::move(roff));
  CHECK_LT(m_result_offsets.size(), std::numeric_limits<roff_index_t>::max())
      << "Need to increase size of roff_index_t";
  return m_result_offsets.size() - 1;
}

void read_cov::cov_pg::add_result(std::unique_ptr<result> r) {
  sub_stats();

  keep_result new_keep;
  new_keep.r = r.get();
  new_keep.keep_until = m_pos + aoffset_t(r->seq_size);

  auto it = std::lower_bound(m_keep_results.begin(), m_keep_results.end(), new_keep);
  if (it != m_keep_results.end()) {
    CHECK(it->r != r.get()) << "Adding an already-present result?";
  }

  m_keep_results.insert(it, std::move(new_keep));

  result_offset roff;
  roff.end_offset = m_pos + aoffset_t(r->seq_size);
  roff.r = r.get();

  roff_index_t idx = add_result_offset(std::move(roff));
  m_unmerged_roffrefs.push_back(idx);

  for (dna_base b : r->a->seq) {
    add_base(b);
    if (m_need_compact) {
      compact_internal();
    }
  }

  flush_pending_reads();

  std::sort(m_paths.begin(), m_paths.end(), paths_comparer());

  flush_old_results();
  m_cov->m_results.emplace(std::move(r));
  add_stats();
}

void read_cov::cov_pg::remove_expired_roffs(path_info& pi, aoffset_t r_size) {
  roff_index_t* new_end =
      std::remove_if(pi.begin(), pi.end(), [this, r_size](roff_index_t roff_idx) {
        DCHECK_LT(roff_idx, m_result_offsets.size());
        const auto& pi_entry = m_result_offsets[roff_idx];

        if (pi_entry.end_offset + r_size <= m_pos) {
          // No more coverage for this one!
          return true;
        }
        return false;
      });
  pi.set_end(new_end);
}

void read_cov::cov_pg::flush_pending_reads() {
  for (auto& prelem : m_pending_reads) {
    const auto& prkey = prelem.first;
    auto& pr = prelem.second;

    auto& roff = m_result_offsets[prkey.roff_idx];
    CHECK(roff.r);
    roff.r->reads.insert(prkey.offset, pr.read_ids, prkey.read_len);
    roff.r->a->read_cov_max_paths = std::max(roff.r->a->read_cov_max_paths, pr.max_paths);
  }
  m_pending_reads.clear();
}

void read_cov::cov_pg::flush_old_results() {
  aoffset_t max_read_len = m_cov->m_options.seqset->max_read_len();

  boost::dynamic_bitset<uint64_t> roff_used(m_result_offsets.size(), false);
  for (auto& path_entry : m_paths) {
    const seqset_range& r = path_entry.first;
    CHECK_LE(r.size(), max_read_len)
        << "Seqset entry size longer than maximum?  This could cause dangling result pointers.";

    auto& pi = path_entry.second;

    for (auto ridx : pi) {
      roff_used[ridx] = true;
    }
  }

  for (roff_index_t roff_idx : m_unmerged_roffrefs) {
    roff_used[roff_idx] = true;
  }

  size_t used_count = roff_used.count();
  // Compact unused roff entries.
  roff_translate_table_t rt;
  rt.resize(roff_used.size());

  std::vector<result_offset> old_result_offsets = std::move(m_result_offsets);
  m_result_offsets.reserve(used_count);

  auto roff_it = old_result_offsets.begin();
  for (size_t i = 0; i != roff_used.size(); ++i) {
    auto& roff = *roff_it;
    if (roff_used[i]) {
      DCHECK(roff.r);
      rt[i] = m_result_offsets.size();
      m_result_offsets.push_back(roff);
    } else {
      roff.r = nullptr;
      rt[i] = std::numeric_limits<roff_index_t>::max();
    }
    ++roff_it;
  }
  CHECK(roff_it == old_result_offsets.end());

  for (auto& p : m_paths) {
    auto& pi = p.second;
    translate_uint32s(pi.begin(), pi.end(), rt);
  }
  translate_uint32s(m_unmerged_roffrefs.data(),
                    m_unmerged_roffrefs.data() + m_unmerged_roffrefs.size(), rt);

  auto keep_results_end =
      std::remove_if(m_keep_results.begin(), m_keep_results.end(), [&](const keep_result& kr) {
        bool expired = (kr.keep_until + max_read_len) < m_pos;
        return expired;
      });
  m_keep_results.erase(keep_results_end, m_keep_results.end());
}

void read_cov::cov_pg::add_base(dna_base b) {
  paths_t old_paths;
  std::swap(old_paths, m_paths);
  m_paths.reserve(old_paths.size());

  // Add the base onto each path
  ++m_pos;
  for (auto& old_path_entry : old_paths) {
    const auto& old_r = old_path_entry.first;
    auto& pi = old_path_entry.second;

    seqset_range new_r = old_r.push_front_drop(b.complement());

    if (new_r.size() != old_r.size() + 1) {
      remove_expired_roffs(pi, new_r.size());
    }

    add_path(new_r, std::move(pi));
  }
  old_paths.clear();

  // Save coverage for any reads we see.  Flat map lets us not
  // allocate a bunch of stuff since we expect to have a reasonably
  // small number of read lengths.
  boost::container::flat_map<int /* read len */, read_id_set> ids_by_lens;

  for (auto& path_entry : m_paths) {
    const auto& r = path_entry.first;
    if (r.size() < m_cov->m_options.readmap->min_read_len()) {
      continue;
    }
    auto read_ids = m_cov->m_options.readmap->entry_to_index_range(r.begin(), r.end());
    ids_by_lens.clear();
    for (uint32_t read_id = read_ids.first; read_id != read_ids.second; ++read_id) {
      readmap::read rd = m_cov->m_options.readmap->get_read_by_id(read_id);
      int read_len = rd.size();

      if (read_len > int(r.size())) {
        continue;
      }
      ids_by_lens[read_len].insert(rd.get_rev_comp().get_read_id());
    }

    auto check_roff_idx = [&](roff_index_t roff_idx, size_t path_priority) {
      auto& roff = m_result_offsets[roff_idx];
      for (const auto& len_and_id : ids_by_lens) {
        int read_len = len_and_id.first;
        const auto& read_ids = len_and_id.second;
        int offset = m_pos - roff.start_offset() - read_len;

        if (offset >= roff.r->seq_size) {
          continue;
        }

        pending_reads_key prkey;
        prkey.roff_idx = roff_idx;
        prkey.offset = offset;
        prkey.read_len = read_len;
        auto& pr = m_pending_reads[prkey];
        pr.read_ids.insert(read_ids);
        pr.max_paths = std::max(pr.max_paths, path_priority);
      }
    };
    const auto& pi = path_entry.second;
    for (roff_index_t roff_idx : pi) {
      check_roff_idx(roff_idx, pi.path_priority());
    }
    for (roff_index_t roff_idx : m_unmerged_roffrefs) {
      check_roff_idx(roff_idx, pi.path_priority());
    }
  }
}

std::unique_ptr<read_cov::cov_pg> read_cov::cov_pg::split() {
  std::unique_ptr<read_cov::cov_pg> new_pg(new read_cov::cov_pg);

  new_pg->m_cov = m_cov;
  new_pg->m_pos = m_pos;

  // When we split we normally immediately add a result, so since
  // we're allocating new storage anyways, include a little extra in
  // the allocation.
  new_pg->m_result_offsets.reserve(m_result_offsets.size() + 1);
  new_pg->m_result_offsets = m_result_offsets;
  CHECK_EQ(new_pg->m_result_offsets.capacity(), new_pg->m_result_offsets.size() + 1);

  new_pg->m_keep_results.reserve(m_keep_results.size() + 1);
  new_pg->m_keep_results = m_keep_results;
  CHECK_EQ(new_pg->m_keep_results.capacity(), new_pg->m_keep_results.size() + 1);

  new_pg->copy_paths_from(count_roffrefs(), m_paths, m_unmerged_roffrefs);

  new_pg->m_max_size = m_max_size;
  new_pg->add_stats();

  return new_pg;
}

void read_cov::cov_pg::adjust_pos(aoffset_t adjust) {
  for (auto& roff : m_result_offsets) {
    roff.end_offset += adjust;
  }
  for (auto& kr : m_keep_results) {
    kr.keep_until += adjust;
  }
  m_pos += adjust;
}

void read_cov::cov_pg::make_roff_translates(roffs_t roffs1, roffs_t roffs2,
                                            roff_translate_table_t& rt1,
                                            roff_translate_table_t& rt2) {
  CHECK(m_result_offsets.empty());
  CHECK(rt1.empty());
  CHECK(rt2.empty());

  rt1.reserve(roffs1.size());
  rt2.reserve(roffs2.size());

  auto in1 = roffs1.begin();
  auto in2 = roffs2.begin();

  while (in1 != roffs1.end() && in2 != roffs2.end()) {
    if (!in1->r) {
      rt1.push_back(std::numeric_limits<roff_index_t>::max());
      ++in1;
      continue;
    }
    if (!in2->r) {
      rt2.push_back(std::numeric_limits<roff_index_t>::max());
      ++in2;
      continue;
    }

    if ((*in1) < (*in2)) {
      auto idx = add_result_offset(std::move(*in1));
      rt1.push_back(idx);
      ++in1;
      continue;
    }

    if ((*in2) < (*in1)) {
      auto idx = add_result_offset(std::move(*in2));
      rt2.push_back(idx);
      ++in2;
      continue;
    }

    CHECK((*in1) == (*in2)) << "In1:\n" << *in1->r->a << "\nin2:\n" << *in2->r->a;

    auto idx = add_result_offset(std::move(*in1));
    rt1.push_back(idx);
    rt2.push_back(idx);
    ++in1;
    ++in2;
  }

  while (in1 != roffs1.end()) {
    if (!in1->r) {
      rt1.push_back(std::numeric_limits<roff_index_t>::max());
      ++in1;
      continue;
    }
    auto idx = add_result_offset(std::move(*in1));
    rt1.push_back(idx);
    ++in1;
  }

  while (in2 != roffs2.end()) {
    if (!in2->r) {
      rt2.push_back(std::numeric_limits<roff_index_t>::max());
      ++in2;
      continue;
    }
    auto idx = add_result_offset(std::move(*in2));
    rt2.push_back(idx);
    ++in2;
  }

  CHECK_EQ(rt1.size(), roffs1.size());
  CHECK_EQ(rt2.size(), roffs2.size());

  CHECK(std::is_sorted(m_result_offsets.begin(), m_result_offsets.end()));
  CHECK(std::adjacent_find(m_result_offsets.begin(), m_result_offsets.end()) ==
        m_result_offsets.end());

  roffs1.clear();
  roffs2.clear();
}

read_cov::cov_pg::path_info read_cov::cov_pg::alloc_big_pi(size_t nelem, size_t prio) {
  m_big_alloc_pi.emplace_back(std::unique_ptr<roff_index_t[]>(new roff_index_t[nelem]));
  m_big_alloc_pi_size += nelem;
  return path_info(m_big_alloc_pi.back().get(), m_big_alloc_pi.back().get() + nelem, prio);
}

read_cov::cov_pg::path_info read_cov::cov_pg::alloc_pi(size_t nelem, size_t prio) {
  if (nelem > k_pi_alloc_chunk) {
    return alloc_big_pi(nelem, prio);
  }

  if (m_alloc_pi_free < nelem || !m_alloc_pi_start) {
    if (!m_cov->m_free_pi_allocs.empty()) {
      m_alloc_pi.emplace_back(std::move(m_cov->m_free_pi_allocs.back()));
      m_cov->m_free_pi_allocs.pop_back();
    } else {
      m_alloc_pi.emplace_back(std::unique_ptr<roff_index_t[]>(new roff_index_t[k_pi_alloc_chunk]));
    }
    m_alloc_pi_start = m_alloc_pi.back().get();
    m_alloc_pi_free = k_pi_alloc_chunk;
  }

  DCHECK_GE(m_alloc_pi_free, nelem);
  DCHECK(m_alloc_pi_start);

  path_info alloced(m_alloc_pi_start, m_alloc_pi_start + nelem, prio);
  m_alloc_pi_start += nelem;
  m_alloc_pi_free -= nelem;
  return alloced;
}

void read_cov::cov_pg::unalloc_pi(roff_index_t* begin, roff_index_t* end) {
  if (m_alloc_pi_start == end) {
    DCHECK_GE(end, begin);
    size_t size_freed = (end - begin);
    m_alloc_pi_free += size_freed;
    m_alloc_pi_start = begin;
  }
}

void read_cov::cov_pg::compact() {
  sub_stats();
  compact_internal();
  add_stats();
}

void read_cov::cov_pg::compact_internal() {
  size_t num_roffrefs = count_roffrefs();

  paths_t old_paths;
  std::swap(old_paths, m_paths);

  std::vector<std::unique_ptr<roff_index_t[]>> old_alloc;
  std::swap(old_alloc, m_alloc_pi);
  m_alloc_pi.clear();
  std::vector<std::unique_ptr<roff_index_t[]>> old_big_alloc;
  std::swap(old_big_alloc, m_big_alloc_pi);
  m_big_alloc_pi.clear();
  m_big_alloc_pi_size = 0;

  m_alloc_pi_free = 0;
  m_alloc_pi_start = nullptr;

  copy_paths_from(num_roffrefs, old_paths, m_unmerged_roffrefs);
  m_unmerged_roffrefs.clear();

  m_cov->m_free_pi_allocs.insert(m_cov->m_free_pi_allocs.end(),
                                 std::make_move_iterator(old_alloc.begin()),
                                 std::make_move_iterator(old_alloc.end()));
  m_need_compact = false;
}

size_t read_cov::cov_pg::count_roffrefs() const {
  size_t tot = 0;
  for (const auto& path_elem : m_paths) {
    const auto& pi = path_elem.second;
    tot += pi.size();
  }

  tot += m_unmerged_roffrefs.size() * m_paths.size();

  return tot;
}

void read_cov::cov_pg::copy_paths_from(size_t num_roffrefs, const paths_t& old_paths,
                                       const std::vector<roff_index_t>& unmerged_roffrefs) {
  CHECK(m_paths.empty());

  path_info compacted = alloc_big_pi(num_roffrefs, path_info::k_prio_not_used);

  auto it = compacted.begin();
  for (const auto& path_elem : old_paths) {
    const auto& pi = path_elem.second;

    auto unmerged_start = std::copy(pi.begin(), pi.end(), it);
    auto unmerged_end =
        std::copy(unmerged_roffrefs.begin(), unmerged_roffrefs.end(), unmerged_start);
    path_info new_pi(it, unmerged_end, pi.path_priority());
    it = unmerged_end;
    m_paths.emplace_back(std::make_pair(path_elem.first, std::move(new_pi)));
  }

  CHECK_EQ(it, compacted.end());
}

void read_cov::cov_pg::join(std::unique_ptr<read_cov::cov_pg> other) {
  sub_stats();
  other->sub_stats();

  CHECK_EQ(m_cov, other->m_cov);

  if (m_pos < other->m_pos) {
    adjust_pos(other->m_pos - m_pos);
  } else if (m_pos > other->m_pos) {
    other->adjust_pos(m_pos - other->m_pos);
  }
  CHECK_EQ(m_pos, other->m_pos);

  paths_t old_paths;
  std::swap(m_paths, old_paths);
  roffs_t old_roffs;
  std::swap(m_result_offsets, old_roffs);

  roff_translate_table_t old_translate;
  roff_translate_table_t new_translate;

  make_roff_translates(std::move(old_roffs), std::move(other->m_result_offsets), old_translate,
                       new_translate);

  translate_uint32s(m_unmerged_roffrefs.data(),
                    m_unmerged_roffrefs.data() + m_unmerged_roffrefs.size(), old_translate);
  translate_uint32s(other->m_unmerged_roffrefs.data(),
                    other->m_unmerged_roffrefs.data() + other->m_unmerged_roffrefs.size(),
                    new_translate);

  m_paths.reserve(old_paths.size() + other->m_paths.size());

  auto it1 = old_paths.begin();
  auto it2 = other->m_paths.begin();

  paths_comparer comparer;
  while (it1 != old_paths.end() && it2 != other->m_paths.end()) {
    if (comparer(*it1, *it2)) {
      add_path(it1->first, std::move(it1->second), &old_translate, &m_unmerged_roffrefs);
      ++it1;
    } else {
      add_path(it2->first, std::move(it2->second), &new_translate, &other->m_unmerged_roffrefs);
      ++it2;
    }
  }

  while (it1 != old_paths.end()) {
    add_path(it1->first, std::move(it1->second), &old_translate, &m_unmerged_roffrefs);
    ++it1;
  }
  while (it2 != other->m_paths.end()) {
    add_path(it2->first, std::move(it2->second), &new_translate, &other->m_unmerged_roffrefs);
    ++it2;
  }

  m_unmerged_roffrefs.clear();
  other->m_unmerged_roffrefs.clear();

  DCHECK(std::is_sorted(m_paths.begin(), m_paths.end(), comparer));

  save_path_priorities();
  if (m_max_size && m_paths.size() > m_max_size) {
    trim_paths_to_max();
    DCHECK(std::is_sorted(m_paths.begin(), m_paths.end(), comparer));
  }

  m_keep_results = join_keep_results(std::move(m_keep_results), std::move(other->m_keep_results));

  other->m_paths.clear();
  other->m_result_offsets.clear();
  other->m_keep_results.clear();

  m_alloc_pi.insert(m_alloc_pi.end(), std::make_move_iterator(other->m_alloc_pi.begin()),
                    std::make_move_iterator(other->m_alloc_pi.end()));
  m_big_alloc_pi.insert(m_big_alloc_pi.end(),
                        std::make_move_iterator(other->m_big_alloc_pi.begin()),
                        std::make_move_iterator(other->m_big_alloc_pi.end()));
  m_big_alloc_pi_size += other->m_big_alloc_pi_size;
  other->m_big_alloc_pi_size = 0;
  other->m_big_alloc_pi.clear();

  other->add_stats();
  add_stats();
}

void read_cov::cov_pg::save_path_priorities() {
  std::vector<size_t> sizes;
  sizes.reserve(m_paths.size());
  for (const auto& p : m_paths) {
    sizes.emplace_back(p.second.size());
  }

  std::sort(sizes.begin(), sizes.end());
  absl::btree_map<size_t /* size */, size_t /* paths cutoff */> cutoffs;
  for (size_t i = 0; i < sizes.size(); ++i) {
    cutoffs[sizes[i]] = i + 1;
  }

  for (auto& p : m_paths) {
    auto it = cutoffs.find(p.second.size());
    CHECK(it != cutoffs.end());
    p.second.increase_path_priority(it->second);
  }
}

void read_cov::cov_pg::trim_paths_to_max() {
  m_cov->on_path_trim(m_paths.size());

  // We've exceeded the maximum number of paths to track.  Trim the
  // paths that have the most result offsets indicating that they
  // were the snarlieiest.

  std::vector<size_t> sizes;
  sizes.reserve(m_paths.size());
  for (const auto& p : m_paths) {
    sizes.emplace_back(p.second.size());
  }

  auto nth = sizes.begin() + m_max_size;
  std::nth_element(sizes.begin(), nth, sizes.end());

  // Trim anything >= this cutoff size.  Don't just cut at the nth
  // element, since it will be arbitrary if there is more than one
  // element with the same size.
  size_t cutoff_size = *nth;
  auto remove_end = std::remove_if(
      m_paths.begin(), m_paths.end(),
      [cutoff_size](const paths_elem_t& elem) { return elem.second.size() >= cutoff_size; });

  m_paths.erase(remove_end, m_paths.end());
}

read_cov::cov_pg::keep_results_t read_cov::cov_pg::join_keep_results(keep_results_t lhs,
                                                                     keep_results_t rhs) {
  auto lhs_it = lhs.begin();
  auto rhs_it = rhs.begin();

  keep_results_t merged;
  merged.reserve(lhs.size() + rhs.size());

  while (lhs_it != lhs.end() && rhs_it != rhs.end()) {
    if (lhs_it->r < rhs_it->r) {
      merged.emplace_back(std::move(*lhs_it));
      ++lhs_it;
      continue;
    }
    if (rhs_it->r < lhs_it->r) {
      merged.emplace_back(std::move(*rhs_it));
      ++rhs_it;
      continue;
    }
    CHECK_EQ(lhs_it->r, rhs_it->r);

    if (lhs_it->keep_until > rhs_it->keep_until) {
      merged.emplace_back(std::move(*lhs_it));
    } else {
      merged.emplace_back(std::move(*rhs_it));
    }
    ++lhs_it;
    ++rhs_it;
  }

  merged.insert(merged.end(), std::make_move_iterator(lhs_it), std::make_move_iterator(lhs.end()));
  merged.insert(merged.end(), std::make_move_iterator(rhs_it), std::make_move_iterator(rhs.end()));

  DCHECK(std::is_sorted(merged.begin(), merged.end()));
  return merged;
}

void read_cov::on_assembly(assembly_ptr a) {
  if (k_cov_debug) {
    std::cout << "Read_cov got assembly: " << *a << "\n";
  }
  if (a->bypass_coverage) {
    if (k_cov_debug) {
      std::cout << "Bypass coverage; skipping\n";
    }
    sort_and_output(std::move(a));
    return;
  }
  if (k_show_stats) {
    ++m_assemblies_seen;
  }
  track_left_offset(a->left_offset);
  advance_to(a->left_offset);
  if (a->left_offset == a->right_offset) {
    m_cur_inserts.emplace_back(std::move(a));
  } else {
    m_cur_non_inserts.emplace_back(std::move(a));
  }
}

std::unique_ptr<read_cov::cov_pg> read_cov::new_pg() {
  std::unique_ptr<cov_pg> pg = make_unique<cov_pg>(this, m_options.seqset->ctx_begin());
  pg->set_max_size(m_options.max_coverage_paths);
  return pg;
}

void read_cov::flush_active_to_here() {
  if (k_cov_debug) {
    std::cout << "Flushing active at " << m_cur_offset << ":  " << m_active.size() << " active, "
              << m_cur_inserts.size() << " cur inserts, " << m_cur_non_inserts.size()
              << " cur non-inserts\n";
  }
  std::unique_ptr<cov_pg> ref_pg;
  std::unique_ptr<cov_pg> rejoin_pg;

  while (!m_active.empty() && m_active.begin()->first == m_cur_offset) {
    std::unique_ptr<cov_pg> pg = std::move(m_active.begin()->second);
    if (rejoin_pg) {
      rejoin_pg->join(std::move(pg));
    } else {
      rejoin_pg = std::move(pg);
    }
    m_active.erase(m_active.begin());
  }
  while (!m_ref_active.empty() && m_ref_active.begin()->first == m_cur_offset) {
    std::unique_ptr<cov_pg> pg = std::move(m_ref_active.begin()->second);
    if (ref_pg) {
      ref_pg->join(std::move(pg));
    } else {
      ref_pg = std::move(pg);
    }
    m_ref_active.erase(m_ref_active.begin());
  }
  if (!m_active.empty()) {
    CHECK_GT(m_active.begin()->first, m_cur_offset);
  }
  if (!m_ref_active.empty()) {
    CHECK_GT(m_ref_active.begin()->first, m_cur_offset);
  }

  if (!ref_pg) {
    ref_pg = new_pg();
  }

  if (!m_cur_inserts.empty()) {
    std::unique_ptr<cov_pg> joined_inserts = ref_pg->split();
    for (auto& a : m_cur_inserts) {
      CHECK_EQ(a->left_offset, a->right_offset);
      auto processed = process_assembly(ref_pg->split(), std::move(a));
      joined_inserts->join(std::move(processed));
    }
    if (rejoin_pg) {
      rejoin_pg->join(std::move(joined_inserts));
    } else {
      rejoin_pg = std::move(joined_inserts);
    }
    m_cur_inserts.clear();
  }

  // Make sure anything we save for later is compacted.
  ref_pg->compact();
  if (rejoin_pg) {
    rejoin_pg->compact();
  }

  for (auto& a : m_cur_non_inserts) {
    if (a->matches_reference) {
      continue;
    }
    aoffset_t right_offset = a->right_offset;
    CHECK_EQ(m_cur_offset, a->left_offset) << *a;
    CHECK_GT(right_offset, m_cur_offset) << *a;

    auto processed = process_assembly(ref_pg->split(), std::move(a));
    auto& active_elem = m_active[right_offset];
    if (active_elem) {
      active_elem->join(std::move(processed));
    } else {
      active_elem = std::move(processed);
      active_elem->compact();
    }
  }

  if (rejoin_pg) {
    rejoin_pg->join(std::move(ref_pg));
  } else {
    rejoin_pg = std::move(ref_pg);
  }
  for (auto& a : m_cur_non_inserts) {
    if (!a) {
      continue;
    }
    CHECK(a->matches_reference);
    aoffset_t right_offset = a->right_offset;
    CHECK_EQ(m_cur_offset, a->left_offset) << *a;
    CHECK_GT(right_offset, m_cur_offset) << *a;

    auto processed = process_assembly(rejoin_pg->split(), std::move(a));
    auto& active_elem = m_ref_active[right_offset];
    if (active_elem) {
      active_elem->join(std::move(processed));
    } else {
      active_elem = std::move(processed);
      active_elem->compact();
    }
  }
  m_cur_non_inserts.clear();
  rejoin_pg.reset();

  size_t max_pi_blocks = k_max_pi_free_blocks_per_active * (m_active.size() + 1);
  if (m_free_pi_allocs.size() > max_pi_blocks) {
    // We probably just finished a difficult merging section or
    // something.  Free up all the extra RAM we had to use
    // temporarily.  Do this before flushing assembly output to other
    // steps which might need to use the RAM.
    m_free_pi_allocs.clear();
  }

  bool do_gc = false;
  if (m_active.empty()) {
    do_gc = true;
  } else {
    time_t now = time(0);
    constexpr int k_gc_interval = 30;
    if (m_last_results_gc + k_gc_interval < now) {
      do_gc = true;
    }
  }

  if (do_gc) {
    m_last_results_gc = time(0);
    gc_results();
    flush_sorted_to(m_cur_offset);
  }
}

void read_cov::gc_results() {
  absl::flat_hash_set<std::unique_ptr<result>> new_results;
  new_results.reserve(m_results.size());
  for (const auto* collection : {&m_active, &m_ref_active}) {
    for (const auto& active_entry : *collection) {
      for (const auto& keep : active_entry.second->keep_results()) {
        auto it = m_results.find(keep.r);
        if (it == m_results.end()) {
          DCHECK(new_results.find(keep.r) != new_results.end());
          continue;
        }
        new_results.insert(m_results.extract(it));
      }
    }
  }
  for (const auto& r : m_results) {
    assembly_ptr a = std::move(r->a);
    aoffset_t left_offset = a->left_offset;
    a->read_coverage.emplace(r->reads.build_and_clear(a->seq.size()));
    sort_and_output(std::move(a));
    untrack_left_offset(left_offset);
  }
  m_results = std::move(new_results);
}

std::unique_ptr<read_cov::cov_pg> read_cov::process_assembly(std::unique_ptr<cov_pg> var_pg,
                                                             assembly_ptr a) {
  std::unique_ptr<result> r = make_unique<result>();
  r->idx = m_result_idx++;
  r->seq_size = a->seq.size();
  r->a = std::move(a);
  r->cov = this;

  var_pg->add_result(std::move(r));

  return var_pg;
}

void read_cov::advance_to(aoffset_t target_offset) {
  while (m_cur_offset < target_offset) {
    flush_active_to_here();

    aoffset_t new_offset = target_offset;

    CHECK_GT(new_offset, m_cur_offset);
    if (!m_active.empty() && new_offset > m_active.begin()->first) {
      new_offset = m_active.begin()->first;
    }
    if (!m_ref_active.empty() && new_offset > m_ref_active.begin()->first) {
      new_offset = m_ref_active.begin()->first;
    }
    CHECK_GT(new_offset, m_cur_offset);

    CHECK(m_cur_inserts.empty());
    CHECK(m_cur_non_inserts.empty());
    m_cur_offset = new_offset;
  }
}

void read_cov::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  flush_sorted();
}

read_cov::read_cov(const assemble_options& options, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_options(options) {}

void read_cov::on_path_trim(size_t paths) {
  if (m_did_notify_path_trim) {
    return;
  }

  m_did_notify_path_trim = true;

  std::stringstream out;

  out << "At position " << m_options.scaffold_name << ":" << m_cur_offset
      << ", read_cov limited trace paths from " << paths << " to " << m_options.max_coverage_paths
      << "; coverage counts may be inaccurate.  Suppressing trim warnings for the rest of "
         "this region.";

  SPLOG_P(LOG_DEBUG, "%s", out.str().c_str());
}

std::ostream& operator<<(std::ostream& os, const read_cov::result_offset& roff) {
  os << "Roff(start=" << roff.start_offset() << ",end=" << roff.end_offset << ")";
  if (roff.r) {
    os << " asm=" << *roff.r->a;
  }
  return os;
}

}  // namespace variants

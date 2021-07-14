#include "modules/build_seqset/part_repo.h"
#include "modules/build_seqset/part_counts.h"
#include "modules/io/make_unique.h"
#include "modules/io/parallel.h"

namespace build_seqset {
namespace {

kmer_t kmer_push_back(kmer_t orig, unsigned kmer_size, dna_base b) {
  return ((orig << 2) | int(b)) & ((1UL << (2 * kmer_size)) - 1);
}

}  // namespace

constexpr unsigned part_repo::k_part_counts_depth;

part_repo::part_repo(unsigned partition_depth,
                     const std::string& ref_path_prefix,
                     const std::string& repo_path)
    : m_depth(partition_depth),
      m_ref_prefix(ref_path_prefix),
      m_repo_path(repo_path) {
  CHECK_GE(m_depth, 1);
  CHECK_LE(m_depth, seq_repository::k_inline_bases);
  flush();
}

dna_sequence part_repo::prefix_for_partition(kmer_t part_num) const {
  return dna_sequence(part_num, m_depth);
}

kmer_t part_repo::partition_for_sequence(const dna_slice& seq) const {
  kmer_t result = 0;
  auto it = seq.begin();
  auto it_end = seq.end();
  for (unsigned i = 0; i < m_depth; ++i) {
    result <<= 2;
    if (it != it_end) {
      result |= int(*it);
      ++it;
    }
  }
  CHECK_LT(result, partition_count());
  return result;
}

void part_repo::add_initial_repo(const dna_slice& reference_data) {
  size_t offset = get_repo_builder()->write_seq(reference_data);
  CHECK_EQ(0, offset);
  flush();
}

void part_repo::write(const dna_slice& seq, unsigned fwd_suffixes,
                      unsigned rc_suffixes) {
  CHECK_GE(seq.size() + 1, fwd_suffixes);
  CHECK_GE(seq.size() + 1, rc_suffixes);

  size_t repo_pos = seq_repository::k_max_offset;
  if (seq.size() > seq_repository::k_inline_bases) {
    repo_pos = get_repo_builder()->write_seq(seq);
  }

  write_using_repo(seq, fwd_suffixes, rc_suffixes, repo_pos);
}

void part_repo::write_using_repo(const dna_slice& seq, unsigned fwd_suffixes,
                                 unsigned rc_suffixes, size_t repo_pos) {
  size_t offset = repo_pos;
  if (repo_pos != seq_repository::k_max_offset) {
    offset = repo_pos + seq_repository::k_inline_bases;
  }
  auto it = seq.begin();
  while (fwd_suffixes) {
    unsigned size = seq.end() - it;
    dna_slice inline_part;
    if (size > seq_repository::k_inline_bases) {
      inline_part = dna_slice(it, seq_repository::k_inline_bases);
    } else {
      inline_part = dna_slice(it, size);
    }

    seq_repository::entry_data e(size, inline_part, offset,
                                 false /* offset is rc */);
    write_raw(e);

    ++it;
    if (size > seq_repository::k_inline_bases) {
      CHECK_LT(offset, seq_repository::k_max_offset);
      ++offset;
    } else {
      offset = seq_repository::k_max_offset;
    }
    fwd_suffixes--;
  }

  it = seq.rcbegin();
  if (seq.size() > seq_repository::k_inline_bases) {
    offset = repo_pos + seq.size() - seq_repository::k_inline_bases;
    CHECK_LT(offset, seq_repository::k_max_offset);
  } else {
    offset = seq_repository::k_max_offset;
  }

  while (rc_suffixes) {
    unsigned size = seq.rcend() - it;
    dna_slice inline_part;
    if (size > seq_repository::k_inline_bases) {
      inline_part = dna_slice(it, seq_repository::k_inline_bases);
    } else {
      inline_part = dna_slice(it, size);
    }

    seq_repository::entry_data e(size, inline_part, offset,
                                 true /* offset is rc */);
    write_raw(e);

    ++it;
    if (size > seq_repository::k_inline_bases) {
      CHECK_GT(offset, 0);
      --offset;
    } else {
      offset = seq_repository::k_max_offset;
    }
    rc_suffixes--;
  }
}

void part_repo::write(const seq_repository::entry_base& e) {
  seq_repository::entry_data d = e.reify_pop();

  write_raw(d);
}

// A thread-local buffer for buffering entry_data writes so we can
// minimize lock contention.
class part_repo::write_entry_buffer : public parallel_local {
 public:
  // Minimum number of entries to accumulate before attempting to flush.
  static constexpr size_t k_flush_size = 1024;
  // Size of buffer before we block on flushing if we can't get a lock.
  static constexpr size_t k_force_flush_size = 2048;

  // If true, output per-thread statistics on lock contention.
  static constexpr bool k_show_flush_sizes = false;
  static constexpr size_t k_num_flush_sizes = (k_force_flush_size / k_flush_size) + 1;

  write_entry_buffer(part_repo* entries) : m_entries(entries) {
    m_bufs.resize(m_entries->m_ref_builders.size());
    if (k_show_flush_sizes) {
      m_flush_size_count.reset(new size_t[k_num_flush_sizes]);
      for (size_t i = 0; i < k_num_flush_sizes; ++i) {
        m_flush_size_count[i] = 0;
      }
    }
  }
  write_entry_buffer() = delete;
  write_entry_buffer(const write_entry_buffer&) = delete;
  write_entry_buffer& operator=(const write_entry_buffer&) = delete;

  void write_entry(const seq_repository::entry_data& data) {
    dna_slice inl = data.inline_bases();

    kmer_t part_num = m_entries->partition_for_sequence(inl);
    auto& buf = m_bufs[part_num];
    if (buf.capacity() < k_flush_size) {
      buf.reserve(k_flush_size);
    }
    buf.emplace_back(data);

    if ((buf.size() % k_flush_size) == 0) {
      bool do_force = false /* default to only flush if we can get the lock */;
      if (buf.size() >= k_force_flush_size) {
        do_force = true;
      }
      flush_part(part_num, do_force);
    }
  }

  void flush() override {
    for (size_t i = 0; i < m_bufs.size(); ++i) {
      flush_part(i, true /* force flush */);
      CHECK(m_bufs[i].empty());
    }

    if (k_show_flush_sizes) {
      std::string out;
      size_t tot = 0;
      for (size_t i = 0; i < k_num_flush_sizes; ++i) {
        tot += m_flush_size_count[i];
      }

      for (size_t i = 0; i < k_num_flush_sizes; ++i) {
        out += " " + printstring("%lu(%.2f%%)", m_flush_size_count[i],
                                 m_flush_size_count[i] * 100. / tot);
      }
      SPLOG("Flush size counts:%s", out.c_str());
    }
  }

  void flush_part(size_t part_num, bool do_force) {
    auto& buf = m_bufs[part_num];
    if (buf.empty()) {
      return;
    }
    size_t orig_size = buf.size();
    m_entries->m_ref_builders[part_num]->write_entries_and_clear(buf, do_force);
    if (do_force) {
      CHECK(buf.empty());
    }
    if (k_show_flush_sizes) {
      if (buf.empty()) {
        m_flush_size_count[orig_size / k_flush_size]++;
      }
    }
  }

 private:
  std::vector<std::vector<seq_repository::entry_data>> m_bufs;
  part_repo* const m_entries;
  std::unique_ptr<size_t[]> m_flush_size_count;
};

constexpr size_t part_repo::write_entry_buffer::k_flush_size;
constexpr size_t part_repo::write_entry_buffer::k_force_flush_size;
constexpr bool part_repo::write_entry_buffer::k_show_flush_sizes;
constexpr size_t part_repo::write_entry_buffer::k_num_flush_sizes;

void part_repo::write_raw(const seq_repository::entry_data& e) {
  auto* st = parallel_pool().get_state();
  if (st) {
    write_entry_buffer* b = st->get_local<write_entry_buffer>(this);
    b->write_entry(e);
  } else {
    CHECK(getenv("TEST_TMPDIR"))
        << "Non-parallel version of part_repo::write_raw should only be called within unit tests.";
    dna_slice inl = e.inline_bases();

    kmer_t part_num = partition_for_sequence(inl);
    m_ref_builders[part_num]->write_entry(e);
  }
}

std::string part_repo::ref_filename(kmer_t part_num,
                                    const std::string& pass_name) const {
  return m_ref_prefix + "-" + pass_name + "-part-" +
         prefix_for_partition(part_num).as_string();
}

void part_repo::for_each_partition(
    const std::string& pass_name,
    const std::function<void(const partition_ref&)>& f,
    progress_handler_t progress) {
  auto parts = partitions(pass_name);
  parallel_for(  //
      0, parts.size(),
      [&](size_t part_num) {
        const partition_ref& part_ref(parts[part_num]);
        f(part_ref);
      },
      subprogress(progress, 0.1, 1));
}

std::vector<part_repo::partition_ref> part_repo::partitions(
    const std::string& pass_name, bool do_pushed, bool delete_on_close) const {
  std::vector<partition_ref> parts;
  parts.resize(partition_count());
  parallel_for(  //
      0, partition_count(), [&](size_t part_num) {
        partition_ref& part_ref(parts[part_num]);

        part_ref.part_id = part_num;
        part_ref.prefix = prefix_for_partition(part_num);
        auto r = open_part_repo(part_num, pass_name);
        r->set_delete_on_close(delete_on_close);
        part_ref.main = std::move(r);
      });

  if (do_pushed) {
    parallel_for(  //
        0, partition_count(), [&](size_t part_num) {
          partition_ref& part_ref(parts[part_num]);
          for (dna_base b : dna_bases()) {
            dna_sequence pushed_prefix;
            pushed_prefix += b;
            pushed_prefix += part_ref.prefix;
            auto& bstart = part_ref.pushed[b].first;
            auto& bend = part_ref.pushed[b].second;

            std::tie(part_ref.pushed_repositories[b], bstart, bend) =
                range_including_prefix(parts, pushed_prefix);
          }
        });
  }
  dna_sequence next_entry;
  unsigned i = partition_count();
  while (i > 0) {
    --i;
    parts[i].next_entry = next_entry;
    if (parts[i].main->begin() != parts[i].main->end()) {
      next_entry = parts[i].main->begin()->sequence();
    }
  }
  return parts;
}

std::tuple<std::shared_ptr<const seq_repository>, seq_repository::iterator,
           seq_repository::iterator>
part_repo::range_including_prefix(const std::vector<partition_ref>& parts,
                                  const dna_slice& seq) const {
  // First, try to get the most specific range we can.
  kmer_t k = 0;
  for (unsigned i = 0; i < m_depth; ++i) {
    if (i < seq.size()) {
      k = kmer_push_back(k, m_depth, seq[i]);
    } else {
      k = kmer_push_back(k, m_depth, dna_base(0));
    }
  }

  auto repo = parts[k].main;
  seq_repository::iterator bstart, bend;

  if (repo->begin() != repo->end()) {
    dna_base last_base = seq[seq.size() - 1];
    if (int(last_base) == 0) {
      bstart = repo->begin();
    } else {
      auto dstart = std::lower_bound(repo->data_begin(), repo->data_end(), seq,
                                     repo->less_than_slice_using_repo());
      bstart = repo->begin() + (dstart - repo->data_begin());
    }

    if (int(last_base) == 3) {
      bend = repo->end();
    } else {
      dna_sequence pushed_end(seq.begin(), seq.end());
      pushed_end[pushed_end.size() - 1] =
          dna_base(1 + int(pushed_end[pushed_end.size() - 1]));
      auto dend =
          std::lower_bound(repo->data_begin(), repo->data_end(), pushed_end,
                           repo->less_than_slice_using_repo());
      bend = repo->begin() + (dend - repo->data_begin());
    }

    if (bstart != bend) {
      return std::make_tuple(repo, bstart, bend);
    }
  }
  return std::make_tuple(nullptr, bstart, bstart);
}

seq_repository::repo_builder* part_repo::get_repo_builder() {
  if (!m_repo_builder) {
    std::lock_guard<std::mutex> l(m_mu);
    if (!m_repo_builder) {
      m_repo_builder.reset(new seq_repository::repo_builder(m_repo_path));
    }
  }
  return m_repo_builder.get();
}

std::shared_ptr<seq_repository> part_repo::open_part_repo(
    kmer_t part_num, const std::string& pass_name) const {
  return make_unique<seq_repository>(ref_filename(part_num, pass_name),
                                     m_repo_slice);
}

void part_repo::dump_part_counts_if_needed() const {
  if (m_part_counts) {
    SPLOG("Pass \"%s\" part counts: %s", m_part_counts_pass_name.c_str(),
          m_part_counts->display_histo().c_str());
  }
}

std::unique_ptr<part_counts> part_repo::release_part_counts(const std::string& pass_name) {
  dump_part_counts_if_needed();
  if (m_part_counts_pass_name == pass_name) {
    return std::move(m_part_counts);
  }
  return nullptr;
}

void part_repo::reset_part_counts(const std::string& pass_name, std::unique_ptr<part_counts> counts) {
  dump_part_counts_if_needed();
  m_part_counts_pass_name = pass_name;
  m_part_counts = std::move(counts);
}

void part_repo::flush() {
  std::lock_guard<std::mutex> l(m_mu);

  bool need_repo_reload = (m_repo_builder || m_repo_slice.size() == 0);

  m_repo_builder.reset();
  m_ref_builders.clear();

  if (need_repo_reload && boost::filesystem::exists(m_repo_path) &&
      boost::filesystem::file_size(m_repo_path) > 0) {
    constexpr bool k_copy_repo_into_ram = true;
    m_repo = membuf();

    if (k_copy_repo_into_ram) {
      mmap_buffer repo_on_disk(m_repo_path);
      m_repo = membuf(new owned_membuf(repo_on_disk.data(), repo_on_disk.size(),
                                       "build_seqset_repo"));
      SPLOG("Loaded %ld bytes of repo into RAM", m_repo.size());
      dna_const_iterator start = dna_const_iterator(
          (const uint8_t*)m_repo.data(), 0 /* offset */, 0 /* not rc */);
      dna_const_iterator end = start + m_repo.size() * 4;
      m_repo_slice = dna_slice(start, end);
    } else {
      m_repo = membuf(new mmap_buffer(m_repo_path));
      dna_const_iterator start = dna_const_iterator(
          (const uint8_t*)m_repo.data(), 0 /* offset */, 0 /* not rc */);
      dna_const_iterator end = start + m_repo.size() * 4;
      m_repo_slice = dna_slice(start, end);
      membuf_cachelist(m_repo).cache_in_memory();
    }
  }
}

void part_repo::open_write_pass(const std::string& pass_name) {
  CHECK(m_ref_builders.empty());
  CHECK_NE(pass_name, "") << "Must specify a write pass name";
  m_ref_builders.resize(partition_count());
  reset_part_counts(pass_name, make_unique<part_counts>(partition_depth() + k_part_counts_depth));
  for (unsigned part_id = 0; part_id < partition_count(); ++part_id) {
    m_ref_builders[part_id] = open_ref_builder(part_id, pass_name);
  }
}

std::unique_ptr<seq_repository::ref_builder> part_repo::open_ref_builder(
    kmer_t part_id, const std::string& pass_name) {
  CHECK_NE(pass_name, "");
  part_counts* counts = nullptr;
  if (pass_name == m_part_counts_pass_name) {
    counts = m_part_counts.get();
  }
  return make_unique<seq_repository::ref_builder>(ref_filename(part_id, pass_name),
                                                  counts);
}

namespace {

struct prefix_searcher {
  bool operator()(const seq_repository::reference& e,
                  const dna_slice& seq) const {
    // Must return true if it is less than the prefix range of "seq".
    dna_compare_result cmp = e.compare_to(seq);
    switch (cmp) {
      case dna_compare_result::FIRST_IS_LESS:
      case dna_compare_result::FIRST_IS_PREFIX:
        return true;
      case dna_compare_result::EQUAL:
      case dna_compare_result::SECOND_IS_LESS:
      case dna_compare_result::SECOND_IS_PREFIX:
        return false;
    }
    return false;
  }

  bool operator()(const dna_slice& seq,
                  const seq_repository::reference& e) const {
    // Must return true if the prefix range of "seq" is less than e.
    dna_compare_result cmp = e.compare_to(seq);
    switch (cmp) {
      case dna_compare_result::SECOND_IS_LESS:
        return true;
      case dna_compare_result::SECOND_IS_PREFIX:
      case dna_compare_result::EQUAL:
      case dna_compare_result::FIRST_IS_LESS:
      case dna_compare_result::FIRST_IS_PREFIX:
        return false;
    }
    return false;
  }
};

}  // namespace

size_t part_repo::write_with_expansions(const seq_repository::entry_base& eb,
                                        unsigned stride, unsigned count) {
  CHECK(count);
  size_t write_count = 1;
  write(eb);
  --count;
  if (!count) {
    return write_count;
  }
  seq_repository::entry e = eb;
  dna_sequence seq_holder = eb.sequence();
  dna_slice seq = seq_holder;
  CHECK(stride);
  unsigned bases_until_output = stride - 1;
  while (seq.size() > 1 && count) {
    seq = seq.subseq(1, seq.size() - 1);

    e = seq_repository::entry(e.pop_front().reify_pop(), m_repo_slice, 0);
    if (bases_until_output) {
      bases_until_output--;
      continue;
    } else {
      count--;
      bases_until_output = stride - 1;
    }
    write(e);
    ++write_count;
  }
  return write_count;
}

}  // namespace build_seqset

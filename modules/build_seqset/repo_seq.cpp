#include "modules/build_seqset/repo_seq.h"
#include "modules/build_seqset/part_counts.h"

namespace build_seqset {

constexpr unsigned seq_repository::k_inline_base_bytes;
constexpr unsigned seq_repository::k_offset_and_rc_bytes;

constexpr unsigned seq_repository::k_inline_bases;
constexpr unsigned seq_repository::k_offset_bits;

constexpr size_t seq_repository::k_max_offset;

seq_repository::seq_repository(const std::string& ref_filename,
                               const std::string& repo_filename)
    : m_ref_filename(ref_filename) {
  if (boost::filesystem::exists(ref_filename) &&
      boost::filesystem::file_size(ref_filename) > 0) {
    m_ref.reset(new mmap_buffer(ref_filename, mmap_buffer::mode::read_write));

    CHECK_EQ(0, m_ref->size() % sizeof(entry_data))
        << "Size " << m_ref->size() << " is not a multiple of "
        << sizeof(entry_data) << " in " << ref_filename;

    m_entry_data_start = reinterpret_cast<entry_data*>(m_ref->mutable_data());
    m_size = m_ref->size() / sizeof(entry_data);
  }

  if (boost::filesystem::exists(repo_filename) &&
      boost::filesystem::file_size(repo_filename) > 0) {
    m_repo.emplace(repo_filename);
    dna_const_iterator start = dna_const_iterator(
        (const uint8_t*)m_repo->data(), 0 /* offset */, 0 /* not rc */);
    dna_const_iterator end = start + m_repo->size() * 4;
    m_repo_slice = dna_slice(start, end);
  }

  SPLOG(
      "Opening sequence repository with %ld entries and a repository of %ld "
      "bases",
      size(), m_repo_slice.size());
}

seq_repository::seq_repository(const std::string& ref_filename,
                               const dna_slice& repo)
    : m_ref_filename(ref_filename) {
  if (boost::filesystem::exists(ref_filename) &&
      boost::filesystem::file_size(ref_filename) > 0) {
    m_ref.reset(new mmap_buffer(ref_filename, mmap_buffer::mode::read_write));

    CHECK_EQ(0, m_ref->size() % sizeof(entry_data))
        << "Size " << m_ref->size() << " is not a multiple of "
        << sizeof(entry_data) << " in " << ref_filename;

    m_entry_data_start = reinterpret_cast<entry_data*>(m_ref->mutable_data());
    m_size = m_ref->size() / sizeof(entry_data);
  }

  m_repo_slice = repo;
}

seq_repository::~seq_repository() {
  if (m_delete_on_close) {
    if (boost::filesystem::exists(m_ref_filename)) {
      boost::filesystem::remove(m_ref_filename);
    }
  }
}

const seq_repository::entry_data* seq_repository::data_begin() const {
  return m_entry_data_start;
}
const seq_repository::entry_data* seq_repository::data_end() const {
  return m_entry_data_start + m_size;
}
seq_repository::entry_data* seq_repository::data_begin() { return m_entry_data_start; }
seq_repository::entry_data* seq_repository::data_end() { return m_entry_data_start + m_size; }

std::function<bool(const seq_repository::entry_data&,
                   const seq_repository::entry_data&)>
seq_repository::less_than_using_repo() const {
  dna_slice seq_repo = repo();
  return [seq_repo](const entry_data& lhs, const entry_data& rhs) -> bool {
    const entry_data& lhs2 = lhs;
    const entry_data& rhs2 = rhs;
    auto cmp = reference(&lhs2, seq_repo)
                   .compare_to(reference(&rhs2, seq_repo));
    switch (cmp) {
      case dna_compare_result::FIRST_IS_LESS:
      case dna_compare_result::FIRST_IS_PREFIX:
        return true;
      default:
        return false;
    }
  };
}

std::function<bool(const seq_repository::entry_data&, const dna_slice&)>
seq_repository::less_than_slice_using_repo() const {
  dna_slice seq_repo = repo();
  return [seq_repo](const entry_data& lhs, const dna_slice& rhs) -> bool {
    switch (reference(&lhs, seq_repo).compare_to(rhs)) {
      case dna_compare_result::FIRST_IS_LESS:
      case dna_compare_result::FIRST_IS_PREFIX:
        return true;
      case dna_compare_result::SECOND_IS_LESS:
      case dna_compare_result::SECOND_IS_PREFIX:
      case dna_compare_result::EQUAL:
        return false;
    }
    return false;
  };
}

seq_repository::iterator seq_repository::begin() const {
  return iterator(m_entry_data_start, repo());
}

seq_repository::iterator seq_repository::end() const {
  return iterator(m_entry_data_start + m_size, repo());
}

void seq_repository::entry_base::check_same_repo(const entry_base& rhs) const {
  CHECK(m_repo.begin() == rhs.m_repo.begin());
  CHECK(m_repo.end() == rhs.m_repo.end());
}

dna_compare_result seq_repository::entry_base::compare_to(const entry_base& rhs) const {
  auto fast_cmp = fast_compare_to(rhs);
  if (fast_cmp) {
    static bool k_double_check_compare = getenv("TEST_TMPDIR") != nullptr;
    if (k_double_check_compare) {
      CHECK_EQ(*fast_cmp, slow_compare_to(rhs));
    }
    return *fast_cmp;
  }
  return slow_compare_to(rhs);
}

boost::optional<dna_compare_result> seq_repository::entry_base::fast_compare_to(
    const entry_base& rhs) const {
  const entry_data* data1 = &get_entry_data();
  unsigned popped1 = popped_count();
  CHECK_LE(popped1, data1->size());
  unsigned size1 = data1->size() - popped1;

  const entry_data* data2 = &rhs.get_entry_data();
  unsigned popped2 = rhs.popped_count();
  CHECK_LE(popped2, data2->size());
  unsigned size2 = data2->size() - popped2;

  if (popped1 < k_inline_bases && popped2 < k_inline_bases) {
    // See if we can order based on just the inline bases so we don't
    // have to look up the full sequence in the repo.
    unsigned compare_len = std::min<unsigned>(
        {k_inline_bases - popped1, size1, k_inline_bases - popped2, size2});
    dna_slice inline1 = data1->inline_bases();
    CHECK_GE(inline1.size(), compare_len);
    inline1 = inline1.subseq(popped1, compare_len);
    dna_slice inline2 = data2->inline_bases();
    CHECK_GE(inline2.size(), compare_len);
    inline2 = inline2.subseq(popped2, compare_len);

    dna_compare_result inline_result = inline1.compare_to(inline2);
    // See if we can order based on just the inline bases.
    switch (inline_result) {
      case dna_compare_result::FIRST_IS_LESS:
      case dna_compare_result::SECOND_IS_LESS:
        return inline_result;

      case dna_compare_result::FIRST_IS_PREFIX:
        if (compare_len == size1) {
          CHECK_GT(size2, compare_len);
          return inline_result;
        }
        break;

      case dna_compare_result::SECOND_IS_PREFIX:
        if (compare_len == size2) {
          CHECK_GT(size1, compare_len);
          return inline_result;
        }
        break;

      case dna_compare_result::EQUAL:
        if (compare_len == size1) {
          if (compare_len == size2) {
            return dna_compare_result::EQUAL;
          } else {
            CHECK_GT(size2, compare_len);
            return dna_compare_result::FIRST_IS_PREFIX;
          }
        } else if (compare_len == size2) {
          CHECK_GT(size1, compare_len);
          return dna_compare_result::SECOND_IS_PREFIX;
        }
        break;
    }
  }

  // Compare repo locations; we could get lucky and not have to
  // dereference them!
  if (data1->size() > k_inline_bases && data2->size() > k_inline_bases) {
    auto repo1_start = get_repo_seq().begin();
    auto repo2_start = rhs.get_repo_seq().begin();
    if (repo1_start.is_rev_comp() == repo2_start.is_rev_comp() &&
        (repo1_start + popped1) == (repo2_start + popped2)) {
      if (size1 < size2) {
        return dna_compare_result::FIRST_IS_PREFIX;
      } else if (size1 > size2) {
        return dna_compare_result::SECOND_IS_PREFIX;
      } else {
        return dna_compare_result::EQUAL;
      }
    }
  }

  return boost::none;
}

dna_compare_result seq_repository::entry_base::slow_compare_to(
    const entry_base& rhs) const {
  check_same_repo(rhs);

  entry_data data1 = reify_pop();
  entry_data data2 = rhs.reify_pop();

  dna_slice inl1 = data1.inline_bases();
  dna_slice inl2 = data2.inline_bases();
  auto inline_cmp = inl1.compare_to(inl2);
  switch (inline_cmp) {
    case dna_compare_result::FIRST_IS_LESS:
    case dna_compare_result::FIRST_IS_PREFIX:
    case dna_compare_result::SECOND_IS_LESS:
    case dna_compare_result::SECOND_IS_PREFIX:
      return inline_cmp;
    case dna_compare_result::EQUAL:
      if (data1.size() <= k_inline_bases) {
        if (data2.size() > k_inline_bases) {
          return dna_compare_result::FIRST_IS_PREFIX;
        } else {
          CHECK_EQ(data1.size(), data2.size());
          return dna_compare_result::EQUAL;
        }
      } else if (data2.size() <= k_inline_bases) {
        CHECK_GT(data1.size(), k_inline_bases);
        return dna_compare_result::SECOND_IS_PREFIX;
      }
      break;
  }
  return reference(&data1, get_repo())
      .get_repo_seq()
      .compare_to(reference(&data2, get_repo()).get_repo_seq());
}

dna_compare_result seq_repository::entry_base::compare_to(
    const dna_slice& rhs) const {
  const entry_data& data = get_entry_data();
  unsigned popped = popped_count();
  CHECK_LE(popped, data.size());

  // See if we can order based on just the inline bases so we don't
  // have to look up the full sequence in the repo.
  if (popped < k_inline_bases) {
    dna_slice inl(dna_const_iterator(data.raw_inline_bases(), 0 /* offset */,
                                     0 /* not rc */) +
                      popped,
                  std::min<unsigned>(data.size(), k_inline_bases) - popped);

    dna_compare_result inline_result = inl.compare_to(rhs);
    switch (inline_result) {
      case dna_compare_result::FIRST_IS_LESS:
      case dna_compare_result::SECOND_IS_LESS:
        return inline_result;

      case dna_compare_result::FIRST_IS_PREFIX:
      case dna_compare_result::EQUAL:
        if (inl.size() == data.size() - popped) {
          return inline_result;
        }
        break;

      case dna_compare_result::SECOND_IS_PREFIX:
        return inline_result;
    }
  }

  // Load from repo.
  return sequence().compare_to(rhs);
}

unsigned seq_repository::entry_base::shared_prefix_length(
    const entry_base& rhs) const {
  check_same_repo(rhs);

  const entry_data* data1 = &get_entry_data();
  unsigned popped1 = popped_count();
  CHECK_LE(popped1, data1->size());
  unsigned size1 = data1->size() - popped1;

  const entry_data* data2 = &rhs.get_entry_data();
  unsigned popped2 = rhs.popped_count();
  CHECK_LE(popped2, data2->size());
  unsigned size2 = data2->size() - popped2;

  if (popped1 < k_inline_bases && popped2 < k_inline_bases) {
    // See if we can order based on just the inline bases so we don't
    // have to look up the full sequence in the repo.
    unsigned compare_len = std::min(std::min(k_inline_bases - popped1, size1),
                                    std::min(k_inline_bases - popped2, size2));

    dna_slice inline1 = data1->inline_bases().subseq(popped1, compare_len);
    dna_slice inline2 = data2->inline_bases().subseq(popped2, compare_len);

    unsigned shared_len = inline1.shared_prefix_length(inline2);
    if ((shared_len < compare_len) || (shared_len == size1) ||
        (shared_len == size2)) {
      return shared_len;
    }
  }

  entry_data reified1;
  if (popped1) {
    reified1 = reify_pop();
    data1 = &reified1;
    popped1 = 0;
  }
  entry_data reified2;
  if (popped2) {
    reified2 = rhs.reify_pop();
    data2 = &reified2;
    popped2 = 0;
  }

  // From this case on, we don't have to worry about popped fronts.
  CHECK_EQ(popped1, 0);
  CHECK_EQ(popped2, 0);

  dna_slice inline1 = data1->inline_bases();
  dna_slice inline2 = data2->inline_bases();

  unsigned inline_shared = inline1.shared_prefix_length(inline2);
  if (inline_shared == size1 || inline_shared == size2) {
    return inline_shared;
  }

  if (inline_shared < k_inline_bases) {
    return inline_shared;
  }

  return inline_shared +
         reference(data1, get_repo())
             .get_repo_seq()
             .shared_prefix_length(reference(data2, get_repo()).get_repo_seq());
}

dna_slice seq_repository::entry_base::get_repo_seq() const {
  const entry_data& data = get_entry_data();
  CHECK_GT(data.size(), k_inline_bases);

  size_t offset;
  bool rc;

  std::tie(offset, rc) = data.offset_and_rc();

  CHECK_LE(offset, m_repo.size());

  dna_const_iterator begin = m_repo.begin() + offset;
  size_t repo_part_size = data.size() - k_inline_bases;
  if (rc) {
    CHECK_GE(offset, repo_part_size);
    begin = begin.rev_comp();
    ++begin;
  } else {
    CHECK_LE(offset + repo_part_size, m_repo.size());
  }

  return dna_slice(begin, repo_part_size);
}

dna_sequence seq_repository::entry_base::sequence() const {
  const entry_data& data = get_entry_data();
  unsigned popped = popped_count();

  dna_slice inl = data.inline_bases();
  unsigned inl_popped = std::min<unsigned>(inl.size(), popped);
  dna_sequence seq(inl.begin() + inl_popped, inl.end());
  popped -= inl_popped;
  if (popped) {
    DCHECK_EQ(0, seq.size());
  }
  if (data.size() > k_inline_bases) {
    dna_slice repo_data = get_repo_seq();
    // TODO(nils): Should there be a dna_sequence::operator+= that takes a
    // dna_slice?
    if (repo_data.size() < popped) {
      popped = repo_data.size();
    }
    seq += dna_sequence(repo_data.begin() + popped, repo_data.end());
  }
  return seq;
}

seq_repository::entry seq_repository::entry_base::pop_front() const {
  const entry_data& data = get_entry_data();
  unsigned popped = popped_count();
  CHECK_GT(data.size(), 0);
  CHECK_LE(popped + 1, data.size());

  return entry(data, m_repo, popped + 1);
}

void seq_repository::entry_data::shift_inline_bases(dna_base new_base) {
  for (unsigned i = 0; i < k_inline_base_bytes - 1; i++) {
    m_inline_bases[i] <<= 2;
    m_inline_bases[i] |= m_inline_bases[i + 1] >> 6;
  }
  m_inline_bases[k_inline_base_bytes - 1] <<= 2;
  m_inline_bases[k_inline_base_bytes - 1] |= int(new_base);

  if (m_size) {
    --m_size;
  }

  if (m_size > k_inline_bases) {
    size_t offset;
    bool rc;
    std::tie(offset, rc) = offset_and_rc();

    if (rc) {
      CHECK_GT(offset, 0);
      --offset;
    } else {
      ++offset;
    }
    set_offset_and_rc(offset, rc);
  }
}

seq_repository::entry_data seq_repository::entry_base::reify_pop() const {
  unsigned popped = popped_count();
  if (!popped) {
    return get_entry_data();
  }

  entry_data d = get_entry_data();

  while (popped) {
    if (d.size() > k_inline_bases) {
      size_t offset;
      bool rc;
      std::tie(offset, rc) = d.offset_and_rc();
      d.shift_inline_bases(reference(&d, get_repo()).get_repo_seq()[0]);
    } else {
      d.shift_inline_bases(dna_base(0));
    }
    --popped;
  }
  return d;
}

seq_repository::entry::entry(const entry_base& orig)
    : entry_base(orig.get_repo()),
      m_data(orig.get_entry_data()),
      m_popped(orig.popped_count()) {}

seq_repository::entry& seq_repository::entry::operator=(
    const entry_base& orig) {
  m_data = orig.get_entry_data();
  m_repo = orig.get_repo();
  m_popped = orig.popped_count();
  return *this;
}

seq_repository::repo_builder::repo_builder(const std::string& filename)
    : m_writer(filename, true /* append */) {
  m_cur_offset = m_writer.pos() * 4;
  m_write_buffer.reset(new unsigned char[k_write_buffer_size]);
  m_write_buffer_it = dna_sequence::iterator(m_write_buffer.get(), 0, false /* not rev comp */);
  m_bases_avail = k_write_buffer_size * 4;
}

seq_repository::repo_builder::~repo_builder() {
  flush_write_buffer(true /* final flush */);
  m_writer.flush();
}

void seq_repository::repo_builder::flush_write_buffer(bool final_flush) {
  CHECK_GE(k_write_buffer_size * 4, m_bases_avail);
  size_t bases_to_write = k_write_buffer_size * 4 - m_bases_avail;
  if (!final_flush) {
    CHECK_EQ(bases_to_write, k_write_buffer_size * 4)
        << "Write buffer should be full before flushing";
  }
  if (!bases_to_write) {
    return;
  }
  size_t bytes_to_write = (bases_to_write + 3) / 4;
  m_writer.write(reinterpret_cast<const char*>(m_write_buffer.get()), bytes_to_write);
  m_write_buffer_it = dna_sequence::iterator(m_write_buffer.get(), 0, false /* not rev comp */);
  m_bases_avail = k_write_buffer_size * 4;
}

void seq_repository::repo_builder::write_seq_unlocked(dna_slice seq) {
  if (!m_bases_avail) {
    flush_write_buffer();
  }

  while (seq.size() > m_bases_avail) {
    m_write_buffer_it = dna_sequence::copy_bases(seq.subseq(0, m_bases_avail), m_write_buffer_it);
    seq = seq.subseq(m_bases_avail, seq.size() - m_bases_avail);
    m_bases_avail = 0;
    flush_write_buffer();
  }

  m_write_buffer_it = dna_sequence::copy_bases(seq, m_write_buffer_it);
  CHECK_GE(m_bases_avail, seq.size());
  m_bases_avail -= seq.size();
}

size_t seq_repository::repo_builder::write_seq(const dna_slice& seq) {
  std::lock_guard<std::mutex> l(m_mu);
  size_t output_offset = m_cur_offset;
  m_cur_offset += seq.size();
  write_seq_unlocked(seq);
  return output_offset;
}

seq_repository::ref_builder::ref_builder(const std::string& filename, part_counts* counts)
    : m_writer(filename, true /* append */), m_part_counts(counts) {
  CHECK_EQ(0, m_writer.pos() % sizeof(entry_data))
      << " Position " << m_writer.pos() << " not a multiple of entry size "
      << sizeof(entry_data);
  m_write_buffer.reserve(k_write_buffer_entries);
}

seq_repository::ref_builder::~ref_builder() { flush(); }

void seq_repository::ref_builder::write_entry(const entry_base& e) {
  write_entry(e.reify_pop());
}

void seq_repository::ref_builder::write_entries_and_clear(std::vector<entry_data>& e,
                                                          bool do_force) {
  if (do_force) {
    {
      std::lock_guard<std::mutex> l(m_mu);
      write_entries_unlocked(e);
    }
    e.clear();
  } else {
    std::unique_lock<std::mutex> l(m_mu, std::try_to_lock);
    if (l.owns_lock()) {
      write_entries_unlocked(e);
      l.unlock();
      e.clear();
    }
  }
}

void seq_repository::ref_builder::write_entries_unlocked(const std::vector<entry_data>& e) {
  if (m_part_counts) {
    for (const auto& data : e) {
      m_part_counts->add(data);
    }
  }
  m_writer.write(reinterpret_cast<const char*>(e.data()), e.size() * sizeof(entry_data));
}

void seq_repository::ref_builder::write_entry(const entry_data& data) {
  std::lock_guard<std::mutex> l(m_mu);

  write_entry_unlocked(data);
}

void seq_repository::ref_builder::write_entry_unlocked(const entry_base& data) {
  write_entry_unlocked(data.reify_pop());
}

void seq_repository::ref_builder::write_entry_unlocked(const entry_data& data) {
  m_write_buffer.push_back(data);
  if (m_write_buffer.size() >= k_write_buffer_entries) {
    flush();
  }
}

void seq_repository::ref_builder::flush() {
  write_entries_unlocked(m_write_buffer);
  m_write_buffer.clear();
}

void seq_repository::ref_builder::write_sequence(
    dna_slice seq, seq_repository::repo_builder& repo, unsigned fwd_suffixes,
    unsigned rc_suffixes) {
  CHECK_GE(seq.size() + 1, fwd_suffixes);
  CHECK_GE(seq.size() + 1, rc_suffixes);

  size_t orig_offset = k_max_offset;
  size_t offset = orig_offset;
  if (seq.size() > k_inline_bases) {
    orig_offset = repo.write_seq(seq);
    offset = orig_offset + k_inline_bases;
  }

  auto it = seq.begin();
  while (fwd_suffixes) {
    unsigned size = seq.end() - it;
    dna_slice inline_part;
    if (size > k_inline_bases) {
      inline_part = dna_slice(it, k_inline_bases);
    } else {
      inline_part = dna_slice(it, size);
    }

    entry_data e(size, inline_part, offset, false /* offset is rc */);
    write_entry(e);

    ++it;
    if (size > k_inline_bases) {
      CHECK_LT(offset, k_max_offset);
      ++offset;
    } else {
      offset = k_max_offset;
    }
    fwd_suffixes--;
  }

  it = seq.rcbegin();
  if (seq.size() > k_inline_bases) {
    offset = orig_offset + seq.size() - k_inline_bases;
    CHECK_LT(offset, k_max_offset);
  } else {
    offset = k_max_offset;
  }

  while (rc_suffixes) {
    unsigned size = seq.rcend() - it;
    dna_slice inline_part;
    if (size > k_inline_bases) {
      inline_part = dna_slice(it, k_inline_bases);
    } else {
      inline_part = dna_slice(it, size);
    }

    entry_data e(size, inline_part, offset, true /* offset is rc */);
    write_entry(e);

    ++it;
    if (size > k_inline_bases) {
      CHECK_LT(offset, k_max_offset);
      CHECK_GT(offset, 0);
      --offset;
    } else {
      offset = k_max_offset;
    }
    rc_suffixes--;
  }
}

bool seq_repository::compare_entry_data_lt(const entry_data& lhs, const entry_data& rhs) const {
  // See if we can compare them without having to resort to the slow lookup of the full sequence.
  int inline_cmp = lhs.inline_bases_cmp(rhs);
  if (inline_cmp || rhs.size() <= k_inline_bases || rhs.size() <= k_inline_bases) {
    return inline_cmp < 0;
  }

  // Inline bases compare as the same; look up full sequence to compare.
  reference lhs_ref(&lhs, repo());
  reference rhs_ref(&rhs, repo());
  dna_compare_result repo_cmp = lhs_ref.get_repo_seq().compare_to(rhs_ref.get_repo_seq());
  switch (repo_cmp) {
    case dna_compare_result::FIRST_IS_LESS:
    case dna_compare_result::FIRST_IS_PREFIX:
      return true;
    default:
      return false;
  }
}

void seq_repository::sort_entry_data(entry_data* begin, entry_data* end) const {
  std::sort(begin, end, [this](const entry_data& lhs, const entry_data& rhs) -> bool {
    return compare_entry_data_lt(lhs, rhs);
  });
}

}  // namespace build_seqset

#include <endian.h>

#include "base/base.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"
#include "modules/io/file_io.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/parallel.h"
#include "modules/io/spiral_file.h"
#include "modules/io/spiral_file_mem.h"

const product_version k_readmap_version{"1.2.0"};
constexpr uint64_t readmap_header::k_magic;
constexpr uint32_t readmap::k_null_index;

readmap::readmap(const std::shared_ptr<seqset>& the_seqset, const std::string& readmap_file_path,
                 const spiral_file_options& sfopts)
    : m_path(readmap_file_path), m_seqset(the_seqset) {
  m_opened.reset(new spiral_file_open_mmap(readmap_file_path, sfopts));
  open_spiral_file(m_opened->open());
}

std::unique_ptr<readmap> readmap::open_anonymous_readmap(const std::string& readmap_file_path) {
  return std::unique_ptr<readmap>(new readmap(readmap_file_path));
}

// Opens an anonymous readmap.  Deprecated.
readmap::readmap(const std::string& readmap_file_path) : m_path(readmap_file_path) {
  m_opened.reset(new spiral_file_open_mmap(readmap_file_path));
  open_spiral_file(m_opened->open());
}

void readmap::open_spiral_file(const spiral_file_open_state& state) {
  state.enforce_max_version("readmap", k_readmap_version);

  m_metadata = state.open_json<readmap_metadata>("readmap.json");
  if (m_seqset) {
    CHECK_EQ(m_seqset->uuid(), m_metadata.seqset_uuid);
  }
  m_sparse_multi.reset(new sparse_multi(state.open_subpart("read_ids")));
  m_read_lengths = int_map_interface::detect_subpart_or_uint8_membuf(state, "read_lengths");

  if (state.subpart_present("mate_loop_ptr")) {
    m_pairing_data_present = true;
    m_mate_loop_ptr = int_map_interface::detect_subpart(state.open_subpart("mate_loop_ptr"));
    m_is_forward.reset(new packed_vector<unsigned, 1>(state.open_subpart("is_forward")));
  } else if (state.subpart_present("mate_pair_ptr")) {
    m_pairing_data_present = true;
    m_mate_pair_ptr.reset(new packed_vector<uint32_t, 32>(state.open_subpart("mate_pair_ptr")));
    m_is_forward.reset(new packed_vector<unsigned, 1>(state.open_subpart("is_forward")));
  }
}

void readmap::enable_mate_loop(
    std::function<dna_sequence(uint64_t /* seqset id */, unsigned /* len */)> lookup_seq,
    progress_handler_t progress) {
  if (m_mate_loop_ptr || !m_mate_pair_ptr) {
    return;
  }
  if (!lookup_seq) {
    lookup_seq = [this](uint64_t seqset_id, unsigned len) {
      return m_seqset->ctx_entry(seqset_id).sequence(len);
    };
  }
  SPLOG("Converting %ld mate pairs to mate loops", m_mate_pair_ptr->size());
  std::unique_ptr<mutable_packed_vector<uint32_t, 32>> mate_loop_ptr(
      new mutable_packed_vector<uint32_t, 32>(m_mate_pair_ptr->size(), "readmap:mate_loops"));

  parallel_for(
      0, mate_loop_ptr->size(),
      [&](size_t start, size_t limit) {
        for (uint32_t read_id = start; read_id != limit; ++read_id) {
          mate_loop_ptr->at(read_id) = std::numeric_limits<uint32_t>::max();
        }
      },
      progress);
  SPLOG("Done initializing mate loops");

  mutable_packed_vector<unsigned, 1> claimed(m_mate_pair_ptr->size(), "mate_pair_to_loop:claimed");

  auto read_id_to_entry = [&](uint32_t read_id, int read_len) -> seqset_range {
    seqset_range r = m_seqset->ctx_entry(index_to_entry(read_id));
    CHECK_LE(read_len, int(r.size()))
        << "Read id: " << read_id << " len: " << read_len << " orig " << get_readlength(read_id)
        << " seq: " << r.begin() << " to " << r.end() << ": " << r.sequence();
    if (int(r.size()) > read_len) {
      r = r.pop_back(r.size() - read_len);
    }
    return r;
  };

  // mate_pair_ptr is a bit messy, in that pairs don't always point
  // back to the original read.  So, we just claim entries in the
  // readmap that have the same length as the read we're converting to
  // a cycle.
  auto claim_read_id = [&](std::pair<uint32_t, uint32_t> range, int read_len, bool forward,
                           bool with_mate) -> uint32_t {
    for (uint32_t read_id = range.first; read_id != range.second; ++read_id) {
      if (has_mate(read_id) != with_mate || get_readlength(read_id) != read_len ||
          get_is_forward(read_id) != forward) {
        continue;
      }
      if (claimed.at(read_id).safe_increment()) {
        continue;
      }

      return read_id;
    }
    LOG(FATAL) << "Unable to claim read id for " << range.first << " to " << range.second;
    return std::numeric_limits<uint32_t>::max();
  };

  parallel_for(
      0, m_seqset->size(),
      [&](size_t start, size_t limit) {
        unsigned dedup_self_pair = 0;
        auto end = m_sparse_multi->iterator_at_source(limit);
        for (auto it = m_sparse_multi->iterator_at_source(start); it != end; ++it) {
          uint64_t seqset_id = (*it).first;
          auto read_id_range = (*it).second;

          for (uint32_t orig_read_id = read_id_range.first; orig_read_id != read_id_range.second;
               ++orig_read_id) {
            if (!get_is_forward(orig_read_id)) {
              continue;
            }
            seqset_range entry = read_id_to_entry(orig_read_id, get_readlength(orig_read_id));
            int read_len = get_readlength(orig_read_id);
            uint32_t orig_mate_read_id = m_mate_pair_ptr->at(orig_read_id);
            if (orig_mate_read_id == readmap::k_null_index) {
              // No mate; 2 element cycle back to self.
              uint32_t read_id =
                  claim_read_id(read_id_range, read_len, true /* forward */, false /* mateless */);
              seqset_range rc_entry =
                  m_seqset->find(lookup_seq(entry.begin(), entry.size()).rev_comp());
              auto rc_read_id_range = entry_to_index_range(rc_entry.begin(), rc_entry.end());
              uint32_t rc_read_id = claim_read_id(rc_read_id_range, read_len, false /* reverse */,
                                                  false /* mateless */);
              CHECK_EQ(mate_loop_ptr->at(read_id), std::numeric_limits<uint32_t>::max());
              mate_loop_ptr->at(read_id) = rc_read_id;
              CHECK_EQ(mate_loop_ptr->at(rc_read_id), std::numeric_limits<uint32_t>::max());
              mate_loop_ptr->at(rc_read_id) = read_id;
              continue;
            }
            uint64_t mate_seqset_id = index_to_entry(orig_mate_read_id);
            int mate_len = get_readlength(orig_mate_read_id);
            // Only process each pair once.
            if (seqset_id < mate_seqset_id) {
              continue;
            }
            if (seqset_id == mate_seqset_id) {
              if (read_len < mate_len) {
                continue;
              }
              if (read_len == mate_len) {
                if (1 & ++dedup_self_pair) {
                  continue;
                }
              }
            }
            uint32_t read_id =
                claim_read_id(read_id_range, read_len, true /* forward */, true /* has mate */);
            uint32_t mate_read_id = claim_read_id(entry_to_index(mate_seqset_id), mate_len,
                                                  true /* forward */, true /* has mate */);
            seqset_range mate_entry = read_id_to_entry(mate_read_id, mate_len);
            seqset_range rc_range =
                m_seqset->find(lookup_seq(entry.begin(), entry.size()).rev_comp());
            uint32_t rc_read_id =
                claim_read_id(entry_to_index_range(rc_range.begin(), rc_range.end()), read_len,
                              false /* reverse */, true /* has mate */);
            CHECK_EQ(get_readlength(rc_read_id), read_len);
            seqset_range rc_mate_range =
                m_seqset->find(lookup_seq(mate_entry.begin(), mate_entry.size()).rev_comp());
            uint32_t rc_mate_read_id =
                claim_read_id(entry_to_index_range(rc_mate_range.begin(), rc_mate_range.end()),
                              mate_len, false /* reverse */, true /* has mate */);
            CHECK_EQ(get_readlength(rc_mate_read_id), mate_len);

            CHECK_EQ(mate_loop_ptr->at(read_id), std::numeric_limits<uint32_t>::max());
            mate_loop_ptr->at(read_id) = rc_read_id;
            CHECK_EQ(mate_loop_ptr->at(rc_read_id), std::numeric_limits<uint32_t>::max());
            mate_loop_ptr->at(rc_read_id) = mate_read_id;
            CHECK_EQ(mate_loop_ptr->at(mate_read_id), std::numeric_limits<uint32_t>::max());
            mate_loop_ptr->at(mate_read_id) = rc_mate_read_id;
            CHECK_EQ(mate_loop_ptr->at(rc_mate_read_id), std::numeric_limits<uint32_t>::max());
            mate_loop_ptr->at(rc_mate_read_id) = read_id;
          }
        }
      },
      progress);
  SPLOG("Done converting mate loops.");
  m_mate_loop_ptr = std::move(mate_loop_ptr);
  m_mate_pair_ptr.reset();
}

bool readmap::get_bit(uint64_t loc) const {
  auto range = m_sparse_multi->lookup(loc);
  return range.first != range.second;
}

int readmap::get_readlength(uint32_t index) const {
  CHECK_LT(index, m_sparse_multi->dest_elem_count());
  return m_read_lengths->get(index);
}

bool readmap::has_mate(uint32_t index) const {
  CHECK(m_pairing_data_present);
  if (m_mate_loop_ptr) {
    uint32_t mate_id = index;
    for (unsigned n = 0; n != 2; ++n) {
      mate_id = m_mate_loop_ptr->get(mate_id);
    }
    return mate_id != index;
  } else {
    return m_mate_pair_ptr->at(index) != std::numeric_limits<uint32_t>::max();
  }
}

uint32_t readmap::get_mate(uint32_t index) const {
  if (!m_pairing_data_present) {
    throw(io_exception("No pairing data present"));
  }
  if (m_mate_loop_ptr) {
    uint32_t rc_read_id = m_mate_loop_ptr->get(index);
    CHECK_NE(rc_read_id, std::numeric_limits<uint32_t>::max());
    uint32_t mate_read_id = m_mate_loop_ptr->get(rc_read_id);
    CHECK_NE(mate_read_id, std::numeric_limits<uint32_t>::max());
    if (mate_read_id == index) {
      throw(io_exception("Read has no mate"));
    }
    return mate_read_id;
  } else {
    if (m_mate_pair_ptr->at(index) == std::numeric_limits<uint32_t>::max()) {
      throw(io_exception("Read has no mate"));
    }
    return m_mate_pair_ptr->at(index);
  }
}

uint64_t readmap::get_mate_entry(uint32_t index) const { return index_to_entry(get_mate(index)); }

bool readmap::get_is_forward(uint32_t index) const {
  CHECK(m_pairing_data_present);
  return m_is_forward->at(index);
}

uint32_t readmap::get_rev_comp(uint32_t index) const {
  if (!m_mate_loop_ptr) {
    throw io_exception{
        "Readmap has no mate loop table; use \"biograph upgrade\" to "
        "construct it"};
  }

  unsigned traverse_count;
  if (get_is_forward(index)) {
    traverse_count = 1;
  } else {
    traverse_count = 3;
  }

  uint32_t read_id = index;
  for (unsigned c = 0; c != traverse_count; ++c) {
    read_id = m_mate_loop_ptr->get(read_id);
    CHECK_NE(read_id, std::numeric_limits<uint32_t>::max());
  }
  return read_id;
}

uint32_t readmap::get_mate_rc(uint32_t index) const {
  if (!m_mate_loop_ptr) {
    throw io_exception{
        "Readmap has no mate loop table; use \"biograph upgrade\" to "
        "construct it"};
  }

  unsigned traverse_count;
  if (get_is_forward(index)) {
    traverse_count = 3;
  } else {
    traverse_count = 1;
  }

  uint32_t read_id = index;
  for (unsigned c = 0; c != traverse_count; ++c) {
    read_id = m_mate_loop_ptr->get(read_id);
    CHECK_NE(read_id, std::numeric_limits<uint32_t>::max());
  }
  return read_id;
}

size_t readmap::get_num_bases() const {
  size_t sum = 0;
  for (size_t i = 0; i < m_read_lengths->size(); i++) {
    sum += m_read_lengths->get(i);
  }
  return sum / 2;
}

readmap::pair_stats readmap::get_pair_stats() const {
  // Calculates and returns number of paired/unpaired reads/bases
  pair_stats result;
  for (size_t i = 0; i < size(); i++) {
    if (has_mate(i)) {
      ++result.paired_reads;
      result.paired_bases += get_readlength(i);
    } else {
      ++result.unpaired_reads;
      result.unpaired_bases += get_readlength(i);
    }
  }

  // We counted each read both forward and backwards, so halve our stats.
  result.paired_reads /= 2;
  result.unpaired_reads /= 2;
  result.paired_bases /= 2;
  result.unpaired_bases /= 2;
  return result;
}

std::vector<int> readmap::fake_coverage(const dna_slice& seq) {
  seqset_range c = m_seqset->ctx_begin();
  std::vector<int> rstart(seq.size());
  std::vector<int> rend(seq.size());
  int pos = 0;
  for (auto it = seq.begin(); it != seq.end(); ++it) {
    dna_base seq_comp = it->complement();
    c = c.push_front_drop(seq_comp);
    if (c.begin() + 1 == c.end() && c.size() == m_seqset->read_len()) {
      int read_len = c.size();
      int start = pos + 1 - read_len;
      if (start < 0) {
        continue;
      }
      rstart[start]++;
      rend[pos]++;
    }
    pos++;
  }
  std::vector<int> r(seq.size());
  int cur = 0;
  for (size_t i = 0; i < seq.size(); i++) {
    cur += rstart[i];
    r[i] = cur;
    cur -= rend[i];
  }
  return r;
}

std::vector<int> readmap::approx_coverage(const dna_slice& seq) const {
  std::vector<std::vector<int>> each = approx_strand_coverage_split(seq);
  for (size_t i = 0; i < seq.size(); i++) {
    each[0][i] += each[1][i];
  }
  return each[0];
}

/*
  Get coverage for a sequence over a specific strand forward or reverse
  const bool forward=True is original 5'->3' sequenced direction
*/
std::vector<int> readmap::approx_strand_coverage(const dna_slice& seq, const bool& forward) const {
  std::vector<std::vector<int>> each = approx_strand_coverage_split(seq);
  int strand = forward ? 0 : 1;
  return each[strand];
}

std::vector<std::vector<int>> readmap::approx_strand_coverage_split(const dna_slice& seq) const {
  seqset_range c = m_seqset->ctx_begin();
  std::vector<std::vector<int>> rstart(2, std::vector<int>(seq.size()));
  std::vector<std::vector<int>> rend(2, std::vector<int>(seq.size()));

  int pos = 0;
  for (auto it = seq.begin(); it != seq.end(); ++it) {
    dna_base seq_comp = it->complement();
    c = c.push_front_drop(seq_comp);
    if (c.begin() + 1 == c.end()) {
      auto loc = c.begin();
      auto range = m_sparse_multi->lookup(loc);
      for (auto index = range.first; index != range.second; index++) {
        int read_len = m_read_lengths->get(index);
        // what we're looking for
        if (read_len > int(c.size())) {
          continue;
        }
        int start = pos + 1 - read_len;
        if (start < 0) {
          continue;
        }
        // We're building the complement, so strand switches here
        int strand = get_is_forward(index) ? 1 : 0;
        rstart[strand][start]++;
        rend[strand][pos]++;
      }
    }
    pos++;
  }

  std::vector<std::vector<int>> ret(2, std::vector<int>(seq.size()));
  int curf = 0;
  int curr = 0;
  for (size_t i = 0; i < seq.size(); i++) {
    curf += rstart[0][i];
    ret[0][i] = curf;
    curf -= rend[0][i];

    curr += rstart[1][i];
    ret[1][i] = curr;
    curr -= rend[1][i];
  }
  return ret;
}

membuf_cachelist readmap::membufs() const {
  membuf_cachelist results;
  if (m_is_forward) {
    results += m_is_forward->membufs();
  }
  if (m_read_lengths) {
    results += m_read_lengths->membufs();
  }
  if (m_mate_loop_ptr) {
    results += m_mate_loop_ptr->membufs();
  }
  if (m_sparse_multi) {
    results += m_sparse_multi->membufs();
  }
  if (m_mate_pair_ptr) {
    results += m_mate_pair_ptr->membufs();
  }
  return results;
}

readmap::read_iterator_range readmap::get_prefix_reads(const seqset_range& r,
                                                       int read_len_limit) const {
  if (m_seqset.get() != r.get_seqset()) {
    throw(io_exception("Cannot use a readmap with a seqset it doesn't belong to."));
  }
  read_len_limit = std::max<int>(read_len_limit, min_read_len());
  if (int(r.size()) < read_len_limit) {
    return read_iterator_range(read_iterator(), read_iterator());
  }
  CHECK_NE(r.begin(), r.end()) << "Invalid seqset range";
  uint32_t initial_read_id = m_sparse_multi->lookup_lower_bound(r.begin());
  return read_iterator_range(
      read_iterator(this, initial_read_id, r.begin(), read_len_limit, r.size()), read_iterator());
}

boost::optional<readmap::read> readmap::get_longest_prefix_read(const seqset_range& r) const {
  boost::optional<uint32_t> read_id = get_longest_prefix_read_id(r);
  if (read_id) {
    return get_read_by_id(*read_id);
  } else {
    return boost::none;
  }
}

boost::optional<uint32_t> readmap::get_longest_prefix_read_id(const seqset_range& r) const {
  CHECK(r.valid());
  // TODO(nils): Make this find reads even if they're before
  // r.begin() or after r.end() like get_prefix_reads does, without
  // being slow.
  boost::optional<uint32_t> result;
  int result_read_len = 0;
  if (r.size() < min_read_len()) {
    return result;
  }
  auto reads = entry_to_index_range(r.begin(), r.end());
  for (uint32_t read_id = reads.first; read_id != reads.second; ++read_id) {
    int read_len = get_readlength(read_id);
    if (read_len > int(r.size())) {
      continue;
    }
    if (read_len > result_read_len) {
      result = read_id;
      result_read_len = read_len;
      if (read_len == int(r.size())) {
        return result;
      }
    }
  }
  return result;
}

readmap::containing_read_iterator_range readmap::get_reads_containing(const seqset_range& r) const {
  if (m_seqset.get() != r.get_seqset()) {
    throw(io_exception("Cannot use a readmap with a seqset it doesn't belong to."));
  }
  if (r.begin() == r.end()) {
    return containing_read_iterator_range(containing_read_iterator(), containing_read_iterator());
  }

  return containing_read_iterator_range(containing_read_iterator(this, r),
                                        containing_read_iterator());
}

readmap::read readmap::get_read_by_id(uint32_t read_id) const {
  if (read_id >= size()) {
    throw(io_exception("Invalid read id"));
  }
  return read(this, read_id);
}

void readmap::calc_read_len_limits() {
  // Limit chunk size so we don't spend all our time processing tiny chunks and locking mutexes.
  constexpr size_t k_min_chunk_size = 1024 * 1024 * 8;
  size_t max_num_chunks = (m_read_lengths->size() / k_min_chunk_size) + 1;

  std::mutex mu;
  auto worklist = make_parallel_for_worklist(
      0, m_read_lengths->size(),
      [&](size_t start, size_t limit) {
        unsigned chunk_min = std::numeric_limits<unsigned>::max();
        unsigned chunk_max = 0;
        for (size_t read_id = start; read_id != limit; ++read_id) {
          unsigned len = m_read_lengths->get(read_id);
          if (len < chunk_min) {
            chunk_min = len;
          }
          if (len > chunk_max) {
            chunk_max = len;
          }
        }
        std::lock_guard<std::mutex> l(mu);
        if (chunk_min < m_min_read_len) {
          m_min_read_len = chunk_min;
        }
        if (chunk_max > m_max_read_len) {
          m_max_read_len = chunk_max;
        }
      },
      max_num_chunks);
  parallel_pool().execute_worklist(worklist);

  m_read_lengths_calculated.store(true, std::memory_order_release);
}

readmap::read_iterator::read_iterator(const readmap* rm, uint32_t read_id, uint64_t seqset_id,
                                      int min_size, int max_size)
    : m_phase(phase::FORWARD),
      m_read(rm, read_id, seqset_id),
      m_min_size(min_size),
      m_max_read_len(max_size),
      m_orig_read_id(read_id),
      m_orig_seqset_id(seqset_id),
      m_orig_max_read_len(max_size) {
  skip_non_matching();
}

void readmap::read_iterator::advance() {
  bool forward = m_phase == phase::FORWARD;
  if (!forward) {
    CHECK_EQ(m_phase, phase::BACKWARD);
  }

  const auto& rm = *m_read.m_readmap;
  if (forward) {
    CHECK_LT(m_read.m_read_id, rm.size());

    ++m_read.m_read_id;
  } else {
    if (m_read.m_read_id == 0) {
      done_direction();
      return;
    }
    --m_read.m_read_id;
  }
}

void readmap::read_iterator::skip_non_matching() {
  while (m_phase != phase::DONE && !skip_non_matching_once()) {
  }
}

bool readmap::read_iterator::skip_non_matching_once() {
  bool forward = m_phase == phase::FORWARD;
  if (!forward) {
    CHECK_EQ(m_phase, phase::BACKWARD);
  }

  const auto& rm = *m_read.m_readmap;
  const auto& ss = *rm.m_seqset;

  bool first_in_group;
  if (forward) {
    if (m_read.m_read_id == rm.size()) {
      done_direction();
      return false;  // More searching needed in other direction
    }

    first_in_group =
        m_read.m_read_id == 0 || rm.m_sparse_multi->dest_is_first_in_group(m_read.m_read_id);
  } else {
    auto prev_read_id = m_read.m_read_id + 1;
    first_in_group =
        prev_read_id == rm.size() || rm.m_sparse_multi->dest_is_first_in_group(prev_read_id);
  }

  if (first_in_group) {
    uint64_t new_seqset_id = rm.m_sparse_multi->reverse_lookup(m_read.m_read_id);
    CHECK_NE(m_read.m_seqset_id, std::numeric_limits<uint64_t>::max());
    while (m_read.m_seqset_id != new_seqset_id) {
      int shared;
      if (forward) {
        CHECK_LT(m_read.m_seqset_id, new_seqset_id) << *this;
        ++m_read.m_seqset_id;
        shared = ss.entry_shared(m_read.m_seqset_id);
      } else {
        CHECK_GT(m_read.m_seqset_id, new_seqset_id) << *this;
        shared = ss.entry_shared(m_read.m_seqset_id);
        --m_read.m_seqset_id;
      }
      CHECK_GT(m_min_size, 0);
      if (shared < int(m_min_size)) {
        done_direction();
        return false;  // More searching needed in other direction
      }
      if (shared < m_max_read_len) {
        m_max_read_len = shared;
      }
    }
  }

  if (rm.get_readlength(m_read.m_read_id) > m_max_read_len) {
    advance();
    return false;  // More searching needed
  }

  // Found a matching read!  No more searching needed.
  return true;
}

void readmap::read_iterator::done_direction() {
  bool forward = m_phase == phase::FORWARD;
  if (!forward) {
    CHECK_EQ(m_phase, phase::BACKWARD);
  }

  if (forward) {
    m_phase = phase::BACKWARD;
    m_read.m_read_id = m_orig_read_id;
    m_read.m_seqset_id = m_orig_seqset_id;
    m_max_read_len = m_orig_max_read_len;

    advance();
  } else {
    m_phase = phase::DONE;
  }
}

std::ostream& operator<<(std::ostream& os, const readmap::read& r) {
  return os << "Read(read_id=" << r.m_read_id << " seqset_id=" << r.get_seqset_id() << ")";
}

std::ostream& operator<<(std::ostream& os, readmap::read_iterator::phase p) {
  using phase = readmap::read_iterator::phase;
  switch (p) {
    case phase::FORWARD:
      return os << "FORWARD";
    case phase::BACKWARD:
      return os << "BACKWARD";
    case phase::DONE:
      return os << "DONE";
    default:
      return os << "INVALID(" << int(p) << ")";
  }
}

std::ostream& operator<<(std::ostream& os, const readmap::read_iterator& it) {
  return os << it.m_phase << ": " << it.m_read << " max=" << it.m_max_read_len
            << " orig=" << it.m_orig_read_id << "," << it.m_orig_seqset_id << ","
            << it.m_orig_max_read_len;
}

readmap::containing_read_iterator::containing_read_iterator(const readmap* rm,
                                                            const seqset_range& r)
    : m_range(r), m_orig_len(r.size()) {
  get_read() = read(rm);
  start_entry();
  skip_non_matching();
}

void readmap::containing_read_iterator::skip_non_matching() {
  while (!at_end()) {
    if (get_read().m_read_id == m_end_read_id) {
      advance_entry();
      continue;
    }

    unsigned readlen = get_read().m_readmap->get_readlength(get_read().m_read_id);
    if (readlen < m_range.size()) {
      advance_read();
      continue;
    }
    return;
  }
}

void readmap::containing_read_iterator::advance_read() {
  CHECK(!at_end());
  CHECK_LT(get_read().m_read_id, m_end_read_id);
  ++get_read().m_read_id;
}

void readmap::containing_read_iterator::start_entry() {
  CHECK(!at_end());
  std::tie(get_read().m_read_id, m_end_read_id) =
      get_read().m_readmap->entry_to_index_range(m_range.begin(), m_range.end());
  get_read().m_seqset_id = m_range.begin();
  m_offset_and_read.first = m_range.size() - m_orig_len;
}

// Depth first search
void readmap::containing_read_iterator::advance_entry() {
  CHECK(!at_end());

  // Try descending
  for (int next_base_int = 0; next_base_int < 4; ++next_base_int) {
    dna_base next_base(next_base_int);

    seqset_range next_pushed = m_range.push_front(next_base);
    if (next_pushed.valid()) {
      m_range = next_pushed;
      start_entry();
      return;
    }
  }

  while (m_range.size() > m_orig_len) {
    // Traverse to next item by incrementing first base
    dna_base prev_base = m_range.front();
    seqset_range popped = m_range.pop_front();
    for (int next_base_int = int(prev_base) + 1; next_base_int < 4; ++next_base_int) {
      dna_base next_base(next_base_int);

      seqset_range next_pushed = popped.push_front(next_base);
      if (next_pushed.valid()) {
        m_range = next_pushed;
        start_entry();
        return;
      }
    }
    m_range = popped;
  }
  m_range = seqset_range();
}

bool readmap::containing_read_iterator::at_end() const { return !m_range.valid(); }

std::ostream& operator<<(std::ostream& os, const readmap::containing_read_iterator& it) {
  os << it.get_read() << " offset=" << it.m_offset_and_read.first << " orig len=" << it.m_orig_len
     << " seqset entry= ";
  if (it.m_range.valid()) {
    os << it.m_range.sequence();
  } else {
    os << " (invalid)";
  }
  return os;
}

uint64_t readmap::mid_to_entry(uint64_t mid_id) const {
  return m_sparse_multi->lookup_mid_to_source(mid_id);
}

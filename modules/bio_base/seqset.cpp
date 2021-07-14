#include "modules/bio_base/seqset.h"
#include "modules/bio_base/readmap.h"
#include "modules/io/version.h"
#include "base/base.h"
#include "modules/io/log.h"
#include "modules/io/parallel.h"
#include "modules/io/msgpack_transfer.h"
#include <algorithm>
#include <random>
#include <queue>

const product_version seqset::seqset_version{"1.1.0"};

size_t seqset::compute_old_format_size(size_t entries) {
  size_t bc_size = bitcount::compute_size(entries);
  return 5 * sizeof(uint64_t) + 4 * bc_size + 2 * sizeof(uint8_t) * entries;
}

seqset::seqset(const spiral_file_create_state &state, size_t entries, unsigned max_entry_len)
    : m_entries(entries), m_is_final(false) {
  state.set_version("seqset", seqset_version);

  seqset_metadata metadata;
  metadata.num_entries = entries;
  state.create_json<seqset_metadata>("seqset.json", metadata);

  m_mutable_fixed = state.create_membuf("fixed", 5 * sizeof(uint64_t));
  m_mutable_entry_sizes = make_unique<mutable_packed_varbit_vector>(
      state.create_subpart("entry_sizes"), entries, max_entry_len);
  m_entry_sizes = m_mutable_entry_sizes->get_int_map_interface();
  m_mutable_shared = make_unique<mutable_packed_varbit_vector>(state.create_subpart("shared"),
                                                               entries, max_entry_len - 1);
  m_shared = m_mutable_shared->get_int_map_interface();

  m_fixed = m_mutable_fixed;

  for (dna_base b : dna_bases()) {
    std::string prev_name = "prev_";
    prev_name += char(b);
    m_prev[b].reset(new bitcount(state.create_subpart(prev_name), entries));
  }

  m_uuid = state.uuid();
}

seqset::seqset(const spiral_file_open_state &state) : m_is_final(true) {
  initialize_from_spiral_file(state);
}

void seqset::initialize_from_spiral_file(const spiral_file_open_state& state) {
  state.enforce_max_version("seqset", seqset_version);

  seqset_metadata metadata = state.open_json<seqset_metadata>("seqset.json");
  m_entries = metadata.num_entries;

  m_fixed = state.open_membuf("fixed");
  CHECK_EQ(5 * sizeof(uint64_t), m_fixed.size());
  m_entry_sizes = int_map_interface::detect_subpart_or_uint8_membuf(state, "entry_sizes");
  CHECK_EQ(m_entries, m_entry_sizes->size());

  m_shared = int_map_interface::detect_subpart_or_uint8_membuf(state, "shared");
  CHECK_EQ(m_entries, m_shared->size());

  for (dna_base b : dna_bases()) {
    std::string prev_name = "prev_";
    prev_name += char(b);
    m_prev[b].reset(new bitcount(state.open_subpart(prev_name)));
  }

  m_uuid = state.uuid();

  compute_read_len();
}

seqset::seqset(const std::string &path, const spiral_file_options &options)
    : m_path(path), m_is_final(true) {
    spiral_file_open_mmap o(path, options);
    initialize_from_spiral_file(o.open());
}

seqset::~seqset() {
  if (m_pop_front_cache) {
    delete[] m_pop_front_cache;
  }
}

void seqset::init() {
  for (dna_base b : dna_bases()) {
    m_prev[b]->init();
  }
}

void seqset::set_shared(size_t row, unsigned shared) {
  CHECK(!m_is_final);
  m_mutable_shared->set(row, shared);
}

void seqset::set_entry_size(size_t row, unsigned new_entry_size) {
  CHECK(!m_is_final);
  m_mutable_entry_sizes->set(row, new_entry_size);
}

void seqset::set_bit(size_t row, dna_base base, bool is_set) {
  CHECK(!m_is_final);
  m_prev[base]->set(row, is_set);
}

bitcount& seqset::mutable_prev(dna_base b) {
  CHECK(!m_is_final);
  return *m_prev[b];
}

void seqset::finalize(progress_handler_t prog) {
  CHECK(!m_is_final);
  size_t offset = 0;
  for (dna_base b : dna_bases()) {
    subprogress(prog, 0.25 * int(b), 0.25 * (int(b) + 1));
    set_fixed(int(b), offset);
    offset += m_prev[b]->finalize(prog);
  }
  set_fixed(dna_base::k_num_bases, offset);
  SPLOG("offset = %zu, m_entries = %zu", offset, m_entries);
  if (offset != m_entries) {
    throw io_exception(printstring("Invalid seqset in finalize: %zu != %zu",
                                   offset, m_entries));
  }
  m_is_final = true;
  compute_read_len();
}

seqset_range seqset::ctx_begin() const {
  seqset_range out(this, 0, 0, 0);
  out.m_end = m_entries;
  return out;
}

class inexact_find_context {
  // Distance is number of mismatches currently
  typedef int dist_t;
  // Location is current range + number of remaining bases
  typedef std::pair<seqset_range, int> location_t;
  // State of search
  dna_slice m_sequence;
};

seqset_range seqset::find(const dna_slice &seq) const {
  seqset_range out(this, 0, 0, m_entries);
  // SPLOG("Seeking: %s", seq.as_string().c_str());
  for (size_t i = 0; i < seq.size(); i++) {
    if (!out.valid())
      break;
    out = out.push_front(seq[seq.size() - 1 - i]);
  }
  return out;
}

seqset_range seqset::find(const dna_sequence &seq) const {
  return find(dna_slice(seq));
}

uint64_t seqset::find_existing(const dna_slice& seq) const {
  uint64_t seqset_id = 0;

  for (size_t i = 0; i < seq.size(); ++i) {
    seqset_id = entry_push_front(seqset_id, seq[seq.size() - i - 1]);
  }

  DCHECK_GE(entry_size(seqset_id), seq.size());

  return seqset_id;
}

uint64_t seqset::find_existing_unique(const dna_slice &seq, size_t expected_unique_len) const {
  while (seq.size() > expected_unique_len) {
    uint64_t seqset_id = find_existing(seq.subseq(0, expected_unique_len));
    uint64_t next_seqset_id = seqset_id + 1;

    // If it's disambiguous, return it.
    if (next_seqset_id == size() || entry_shared(next_seqset_id) < expected_unique_len) {
      return seqset_id;
    }

    // Otherwise, try doubling the unique len.
    expected_unique_len *= 2;
  }

  return find_existing(seq);
}

static bool find_near_recursive(std::vector<seqset_range> &out,
                                const dna_slice &seq, size_t max_mismatch,
                                size_t max_results, const seqset_range &cur,
                                int seq_pos) {
  /*
  SPLOG("Cur = [%lu,%lu), pos = %d, max_mismatch = %zu",
          cur.begin(), cur.end(), seq_pos, max_mismatch);
  */
  // Check if I made it though whole sequence
  if (seq_pos == -1) {
    // Room to push?
    if (out.size() >= max_results) {
      return false; // Nope, return failure
    }
    // Yup, add to the list
    out.push_back(cur);
    return true;
  }
  // Simple case for no mismatches left
  if (max_mismatch == 0) {
    seqset_range next = cur.push_front(seq[seq_pos]);
    if (!next.valid()) {
      return true;
    }
    return find_near_recursive(out, seq, 0, max_results, next, seq_pos - 1);
  }
  // Try all the cases
  for (dna_base b : dna_bases()) {
    seqset_range next = cur.push_front(b);
    if (!next.valid()) {
      continue;
    }
    int new_mismatch = max_mismatch - ((b == seq[seq_pos]) ? 0 : 1);
    bool r = find_near_recursive(out, seq, new_mismatch, max_results, next,
                                 seq_pos - 1);
    if (!r) {
      return r;
    }
  }
  return true;
}

bool seqset::find_near(std::vector<seqset_range> &out, const dna_slice &seq,
                       size_t max_mismatch, size_t max_results) const {
  if (max_mismatch == 0) {
    // Special case 0 mismatches
    seqset_range r = find(seq);
    if (max_results == 0) {
      return false;
    }
    if (r.valid()) {
      out.push_back(r);
    }
    return true;
  }
  return find_near_recursive(out, seq, max_mismatch, max_results, ctx_begin(),
                             seq.size() - 1);
}

dna_base seqset::entry_get_base(uint64_t offset) const {
  // unrolled 'binary search'
  return dna_base((offset < get_fixed(2))
                      ? (offset < get_fixed(1) ? 0 : 1)
                      : (offset < get_fixed(3) ? 2 : 3));
}

void seqset::compute_read_len() {
  seqset_range c = ctx_entry(0);
  bool good = true;
  while (good) {
    good = false;
    for (dna_base b : dna_bases()) {
      seqset_range n = c.push_front(b);
      if (!n.valid()) {
        continue;
      }
      good = true;
      c = n;
      break;
    }
  }
  m_read_len = c.size();
}

void seqset::populate_pop_front_cache(progress_handler_t progress) const {
  CHECK(!m_pop_front_cache);
  SPLOG("seqset::populate_pop_front_cache>, m_entries = %zu", m_entries);

  SPLOG("seqset::populate_pop_front_cache> doing resize");
  CHECK_LT(m_entries, (uint64_t(std::numeric_limits<uint8_t>::max())
                       << 32) + std::numeric_limits<uint32_t>::max())
      << "Too many entries to fit in 5 bytes of pop front cache";
  m_pop_front_cache = new uint8_t[m_entries * 5];

  SPLOG("seqset::populate_pop_front_cache> resize done");
  parallel_for(  //
      0, m_entries,
      [this](size_t start, size_t limit) {
        dna_base_array<uint64_t> base_offset;
        for (dna_base b : dna_bases()) {
          base_offset[b] = get_fixed(int(b)) + m_prev[b]->count(start);
        }

        for (uint64_t i = start; i < limit; i++) {
          for (dna_base b : dna_bases()) {
            if (m_prev[b]->get(i)) {
              size_t offset = base_offset[b]++;

              uint8_t* entry = &m_pop_front_cache[offset * 5];
              *entry = i >> 32;
              ++entry;
              *reinterpret_cast<uint32_t*>(entry) = i;
            }
          }
        }
        for (dna_base b : dna_bases()) {
          CHECK_EQ(get_fixed(int(b)) + m_prev[b]->count(limit), base_offset[b]);
        }
      },
      progress);
  SPLOG("seqset::populate_pop_front_cache> population complete");
}

seqset_range seqset_range::next() const {
  seqset_range out(m_seqset, m_seq_size, m_end, m_end);
  out.inner_next();
  return out;
}

seqset_range seqset_range::push_front(const dna_base &b) const {
  // Make sure all is valid
  if (!valid()) {
    throw io_exception("Cannot push_front on to an invalid k-mer");
  }
  // Subtend begin and end within base
  uint64_t sub_begin = m_seqset->m_prev[b]->count(m_begin); // count is O(a,i)
  uint64_t sub_end = m_seqset->m_prev[b]->count(m_end);     // count is O(a,i)
  // Pick out fixed component
  uint64_t fixed =
      m_seqset->get_fixed(int(b)); // m_fixed is C(a), beginning of range for a
  // Offset by fixed component
  uint64_t new_begin = fixed + sub_begin;
  uint64_t new_end = fixed + sub_end;
  // Check if we need to kick begin forward 1
  if (new_begin < new_end && m_seqset->entry_size(new_begin) < m_seq_size + 1)
    new_begin++;
  // Return results
  return seqset_range(m_seqset, m_seq_size + 1, new_begin, new_end);
}

void seqset::init_shared_lt_search() {
  static std::mutex g_mu;
  std::lock_guard<std::mutex> l(g_mu);
  if (m_shared_lt_search) {
    return;
  }

  m_shared_lt_search = make_unique<less_than_search>(m_shared.get());
}

seqset_range seqset_range::push_front_drop(const dna_base &b,
                                           unsigned min_ctx) const {
  // Make sure all is valid
  CHECK(valid());
  // Get base as array index
  int i = (int)b;
  // Pick out fixed component
  uint64_t fixed = m_seqset->get_fixed(i);
  // Subtend begin and end within base
  uint64_t o_begin = m_begin;
  uint64_t o_end = m_end;
  unsigned o_context = m_seq_size;
  uint64_t sub_begin = m_seqset->m_prev[b]->count(o_begin);
  uint64_t sub_end = m_seqset->m_prev[b]->count(o_end);

  if (o_context < min_ctx) {
    return seqset_range(m_seqset, 0, 0, 0);
  }
  while (sub_begin == sub_end ||
         (sub_begin + 1 == sub_end &&
          m_seqset->entry_size(fixed + sub_begin) < o_context + 1)) {
    // SPLOG("sub_begin = %lu, sub_end = %lu, bpos = %d, o_context = %d",
    //	sub_begin, sub_end, m_seqset->entry_sizes(fixed + sub_begin),
    // o_context);
    unsigned drop = std::max(m_seqset->entry_shared(o_begin),
                            (o_end == m_seqset->m_entries
                                 ? unsigned(0)
                                 : m_seqset->entry_shared(o_end)));
    if (sub_begin != sub_end) {
      drop =
          std::max(drop, unsigned(m_seqset->entry_size(fixed + sub_begin) - 1));
    }
    if (drop < min_ctx) {
      return seqset_range(m_seqset, 0, 0, 0);
    }
    bool update_begin = false;
    bool update_end = false;

    if (!m_seqset->m_shared_lt_search) {
      const_cast<seqset *>(m_seqset)->init_shared_lt_search();
    }

    if (o_begin > 0 && m_seqset->entry_shared(o_begin) >= drop) {
      uint64_t drop_begin = m_seqset->m_shared_lt_search->next_backward_lt(o_begin, drop);
      CHECK_LT(drop_begin, o_begin);
      o_begin = drop_begin;
      if (o_begin > 0) {
        DCHECK_LT(m_seqset->entry_shared(o_begin), drop);
      }
      update_begin = true;
    }
    if (o_end < m_seqset->m_entries &&
           m_seqset->entry_shared(o_end) >= drop) {
      uint64_t drop_end = m_seqset->m_shared_lt_search->next_forward_lt(o_end, drop);
      CHECK_GT(drop_end, o_end);
      o_end = drop_end;
      if (o_end < m_seqset->m_entries) {
        DCHECK_LT(m_seqset->entry_shared(o_end), drop);
      } else {
        DCHECK_EQ(o_end, m_seqset->m_entries);
      }
      update_end = true;
    }
    if (update_begin) {
      sub_begin = m_seqset->m_prev[b]->count(o_begin);
    }
    if (update_end) {
      sub_end = m_seqset->m_prev[b]->count(o_end);
    }
    o_context = drop;
  }
  // Offset by fixed component
  uint64_t new_begin = fixed + sub_begin;
  uint64_t new_end = fixed + sub_end;
  // Check if we need to kick begin forward 1
  if (new_begin < new_end && m_seqset->entry_size(new_begin) < o_context + 1)
    new_begin++;
  // Return results
  return seqset_range(m_seqset, o_context + 1, new_begin, new_end);
}

bool seqset_range::find_maximal_prefix_reads(std::set<seqset_range> &results, uint32_t max_reads,
								 unsigned min_overlap, const seqset_bitmap_base &read_bitmap) const
{
	if (read_bitmap.get_bit(m_begin) && is_maximal()) {
		if (results.size() < max_reads) {
			//results.push_back(*this);
			results.insert(*this);
		} else {
			return false;
		}
	}

	for (dna_base b : dna_bases()) {
		seqset_range added_prefix = push_front_drop(b, min_overlap);
		if (!added_prefix.valid()) {
			continue;
		}
		if (!added_prefix.find_maximal_prefix_reads(results, max_reads, min_overlap + 1, read_bitmap)) {
      		return false;
		}
	}

	return true;
}

bool seqset_range::is_maximal() const {
  if (m_begin + 1 != m_end) {
    return false;
  }
  if (size() != m_seqset->entry_size(m_begin)) {
    return false;
  }
  for (dna_base b : dna_bases()) {
    seqset_range pushed_range = push_front(b);
    if (pushed_range.valid()) {
      return false;
    }
  }
  return true;
}

bool seqset_range::find_full_prefix_reads(std::vector<seqset_range> &results, uint32_t max_reads,
								 unsigned min_overlap, const readmap &read_bitmap) const
{
	if (read_bitmap.get_bit(m_begin) && is_full_read(read_bitmap)) {
		if (results.size() < max_reads) {
			results.push_back(*this);
			//results.insert(*this);
		} else {
			return false;
		}
	}

	for (dna_base b : dna_bases()) {
		seqset_range added_prefix = push_front_drop(b, min_overlap);
		if (!added_prefix.valid()) {
			continue;
		}
		if (!added_prefix.find_full_prefix_reads(results, max_reads, min_overlap + 1, read_bitmap)) {
      		return false;
		}
	}

	return true;
}

bool seqset_range::is_full_read(const readmap &read_bitmap) const {
  auto index_range = read_bitmap.entry_to_index(m_begin);
  for (auto idx = index_range.first; idx != index_range.second; idx++) {
    if(read_bitmap.get_readlength(idx) == int(size())) {
      return true;
    }
  }
  return false;
}


bool seqset_range::find_overlap_reads(overlaps_t &results, uint32_t max_reads,
                                      unsigned min_overlap,
                                      const seqset_bitmap_base &read_bitmap,
                                      bool rely_on_read_bitmap,
                                      unsigned added) const {
  if (added != 0 &&
      (rely_on_read_bitmap
           ? read_bitmap.get_bit(m_begin) && m_begin + 1 == m_end &&
                 size() == m_seqset->entry_size(m_begin)
           : read_bitmap.get_bit(m_begin) && is_maximal())) {
    if (results.size() < max_reads) {
      results.emplace(m_begin, size() - added);
      return true;
    }
    return false;
  }

  for (dna_base b : dna_bases()) {
    seqset_range added_prefix = push_front_drop(b, min_overlap);
    if (!added_prefix.valid()) {
      continue;
    }
    if (!added_prefix.find_overlap_reads(results, max_reads, min_overlap + 1,
                                         read_bitmap, rely_on_read_bitmap, added + 1)) {
      return false;
    }
  }

  return true;
}

seqset_range seqset::read_ctx_entry(const readmap& rm, uint32_t readentry) const {
  auto idx = rm.index_to_entry(readentry);
  auto seq = ctx_entry(idx);
  return seq.pop_back(seq.size() - rm.get_readlength(readentry));
}

namespace {

struct overlap_queue_entry {
  seqset_range range;
  unsigned overlap_bases;
  unsigned added;

  bool operator<(const overlap_queue_entry &rhs) const {
    return overlap_bases < rhs.overlap_bases;
  }
};

}  // namespace

std::vector<overlap_result_t> seqset_range::find_overlap_reads_fair(
    uint32_t max_overlaps, unsigned min_overlap,
    const seqset_bitmap_base &read_bitmap, bool rely_on_read_bitmap,
    unsigned added) const {
  std::vector<overlap_result_t> results;
  std::priority_queue<overlap_queue_entry> queue;
  queue.push({.range = *this,
          .overlap_bases = size(),
          .added = 0});

  while (!queue.empty()) {
    overlap_queue_entry entry = queue.top();
    queue.pop();

    if (entry.overlap_bases < min_overlap) {
      return results;
    }

    if (entry.added &&
        (rely_on_read_bitmap
             ? read_bitmap.get_bit(entry.range.begin()) &&
                   entry.range.begin() + 1 == entry.range.end() &&
                   entry.range.size() ==
                       m_seqset->entry_size(entry.range.begin())
             : read_bitmap.get_bit(entry.range.begin()) &&
                   entry.range.is_maximal())) {
      overlap_result_t overlap;
      overlap.seqset_id = entry.range.begin();
      overlap.overlap_bases = entry.overlap_bases;
      results.push_back(overlap);
      if (results.size() > max_overlaps) {
        return results;
      }
      continue;
    }

    for (dna_base b : dna_bases()) {
      overlap_queue_entry new_entry;
      new_entry.range = entry.range.push_front_drop(b, entry.added + min_overlap);
      if (!new_entry.range.valid()) continue;
      new_entry.added = entry.added + 1;
      DCHECK_GT(new_entry.range.size(), new_entry.added);
      new_entry.overlap_bases = new_entry.range.size() - new_entry.added;
      DCHECK_GE(new_entry.overlap_bases, min_overlap);
      queue.push(new_entry);
    }
  }
  return results;
}

seqset_range seqset_range::pop_front() const {
  if (!valid()) {
    throw io_exception("Cannot pop_front from an invalid k-mer");
  }
  if (m_seq_size == 0) {
    throw io_exception("Cannot pop_front from an empty k-mer");
  }
  dna_base b = front();
  unsigned new_context = m_seq_size - 1;
  uint64_t new_begin = inner_pop_front(b, m_begin);
  uint64_t new_end = new_begin + 1;
  while (new_begin && m_seqset->entry_shared(new_begin) >= new_context) {
    new_begin--;
  }
  while (new_end < m_seqset->m_entries &&
         m_seqset->entry_shared(new_end) >= new_context) {
    new_end++;
  }
  return seqset_range(m_seqset, new_context, new_begin, new_end);
}

seqset_range seqset_range::pop_back(size_t count) const {
  if (!valid()) {
    throw io_exception("Cannot pop_back from an invalid k-mer");
  }
  if (m_seq_size < count) {
    throw io_exception(
        boost::format("Cannot pop_back %1% bases from an k-mer of length %2%") %
        count % static_cast<int>(m_seq_size));
  }
  unsigned new_context = m_seq_size - count;
  uint64_t new_begin = m_begin;
  uint64_t new_end = m_end;
  // TODO: This could be done via log(n) method via a 'min_tree'
  while (new_begin && m_seqset->entry_shared(new_begin) >= new_context)
    new_begin--;
  while (new_end < m_seqset->m_entries &&
         m_seqset->entry_shared(new_end) >= new_context)
    new_end++;
  return seqset_range(m_seqset, new_context, new_begin, new_end);
}

seqset_range seqset_range::truncate(size_t count) const {
  if (!valid()) {
    throw io_exception("Cannot truncate from an invalid k-mer");
  }

  if (size() > count) {
    return pop_back(size() - count);
  } else {
    return *this;
  }
}

dna_base seqset_range::front() const {
  if (!valid()) {
    throw io_exception("Cannot call front on an invalid k-mer");
  }
  if (m_seq_size == 0) {
    throw io_exception("Cannot call front on an empty k-mer");
  }
  return m_seqset->entry_get_base(m_begin);
}

dna_sequence seqset_range::sequence(int size) const {
  if (!valid()) {
    throw io_exception("Cannot call sequence on an invalid k-mer");
  }
  dna_sequence tot;
  uint64_t cur = m_begin;
  if (size < 0 || size > int(m_seq_size)) {
    size = m_seq_size;
  }
  tot.reserve(size);
  for (int i = 0; i < size; i++) {
    dna_base b = m_seqset->entry_get_base(cur);
    tot.push_back(b);
    cur = inner_pop_front(b, cur);
  }
  return tot;
}

void seqset_range::inner_next() {
  // Find first entry >= begin with sufficent context
  while (m_seqset->entry_size(m_begin) < m_seq_size &&
         m_begin < m_seqset->m_entries) {
    m_begin++;
  }
  // If I hit the end, make invalid, otherwise, set end to next element
  if (m_begin == m_seqset->m_entries) {
    m_end = m_begin;
  } else {
    m_end = m_begin + 1;
  }
  // While end overlaps with begin, move end forward
  while (m_seqset->entry_shared(m_end) >= m_seq_size &&
         m_end < m_seqset->m_entries) {
    m_end++;
  }
}

uint64_t seqset_range::inner_pop_front(dna_base b, uint64_t offset) const {
  if (m_seqset->is_pop_front_cached()) {
    return m_seqset->entry_pop_front(offset);
  }

  uint64_t base_offset = offset - m_seqset->get_fixed(int(b)) + 1;
  auto itBegin = m_seqset->m_prev[b]->begin();
  auto itEnd = m_seqset->m_prev[b]->end();
  auto it = std::lower_bound(itBegin, itEnd, base_offset);
  uint64_t tot_offset = it - itBegin;
  return tot_offset - 1;
}

unsigned seqset_range::shared_prefix_length(const seqset_range &rhs) const {
  CHECK(valid());
  CHECK(rhs.valid());
  unsigned shared_bases = std::min<unsigned>(size(), rhs.size());
  if (end() > rhs.begin() && rhs.end() > begin()) {
    // these overlap, so one is the prefix of the other.
    return shared_bases;
  }

  uint64_t shared_start, shared_end;
  if (end() > rhs.begin()) {
    shared_start = rhs.end();
    shared_end = begin();
  } else {
    shared_start = end();
    shared_end = rhs.begin();
  }

  for (uint64_t idx = shared_start; idx <= shared_end; ++idx) {
    unsigned s = m_seqset->entry_shared(idx);
    if (s < shared_bases) {
      shared_bases = s;
    }
  }
  return shared_bases;
}

boost::optional<uint64_t> seqset::find_kmer(const dna_slice &seq) const {
  CHECK_GT(seq.size(), 1);

  seqset_range range = find(seq);
  if (!range.valid()) {
    return boost::none;
  }

  return range.begin();
}

boost::optional<uint64_t> seqset::kmer_push_front(uint64_t seqset_id,
                                                  unsigned kmer_size,
                                                  dna_base base) const {
  DCHECK_LT(entry_shared(seqset_id), kmer_size);
  DCHECK_GE(entry_size(seqset_id), kmer_size);

  // We want to make sure there exists a front associated with this
  // kmer without the last base.  First, try linear searching
  // backward.
  uint64_t backward = seqset_id;
  while (backward && !entry_has_front(backward, base) &&
         entry_shared(backward) >= (kmer_size - 1)) {
    --backward;
  }
  if (entry_has_front(backward, base)) {
    seqset_id = backward;
  } else {
    // Otherwise, try linear searching forward; if we don't find a front, there
    // isn't one.
    uint64_t forward = seqset_id;
    while (!entry_has_front(forward, base)) {
      ++forward;
      if (forward == size() || entry_shared(forward) < (kmer_size - 1)) {
        return boost::none;
      }
    }
    seqset_id = forward;
  }

  DCHECK(entry_has_front(seqset_id, base));
  uint64_t pushed = entry_push_front(seqset_id, base);

  if (entry_size(pushed) < kmer_size) {
    // There is a front, but it doesn't overlap enough to use.
    return boost::none;
  }

  // Since we may not be on the first entry that starts with the new pushed
  // kmer, linear search back until we find it.
  uint64_t walk_back = pushed;
  while (walk_back && entry_shared(walk_back) >= kmer_size) {
    --walk_back;
  }

  uint64_t new_entry = walk_back;

  DCHECK_GE(entry_size(new_entry), kmer_size);
  DCHECK_LT(entry_shared(new_entry), kmer_size);
  return new_entry;
}

membuf_cachelist seqset::membufs() const {
  membuf_cachelist results{m_shared->membufs(), m_entry_sizes->membufs()};
  for (dna_base b : dna_bases()) {
    results += m_prev[b]->membufs();
  }
  return results;
}

namespace {

// Header for old non-zip file format.  Deprecated.
struct seqset_header {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(version);
    FIELD(ref_size);  // TODO: Validate ref matched
    FIELD(seqset_offset);
    FIELD(refbits_offset);
  };
  product_version version;
  uint64_t ref_size;  // If 0, no reference data
  uint64_t seqset_offset;
  uint64_t refbits_offset;
};

}  // namespace

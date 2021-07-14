#include "modules/bio_base/kmer.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/mapred/resource_manager.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/log.h"
#include "modules/io/parallel.h"
#include "modules/io/track_mem.h"

class kmer_set::kmer_tail_reference
    : public boost::totally_ordered<
          boost::less_than_comparable<kmer_tail_reference, kmer_t>> {
 public:
  kmer_tail_reference() = delete;
  kmer_tail_reference(const kmer_tail_reference&) = default;
  kmer_tail_reference(kmer_tail_reference&&) = default;
  kmer_tail_reference(kmer_set* ks, lookup_t index)
      : m_ks(ks), m_table_index(index) {
  }

  operator kmer_t() const {
    kmer_t ret = 0;
    unsigned char* pos = get_pos();
    for (unsigned i = 0; i < m_ks->m_tail_bytes; i++) {
      ret <<= 8;
      ret |= kmer_t(pos[i]);
    }
    return ret;
  }

  kmer_tail_reference& operator=(const kmer_tail_reference& rhs) {
    DCHECK(m_ks);
    DCHECK_EQ(m_ks, rhs.m_ks);
    unsigned char* pos = get_pos();
    unsigned char* rhs_pos = rhs.get_pos();
    if (pos != rhs_pos) {
      memcpy(pos, rhs_pos, m_ks->m_tail_bytes);
    }
    return *this;
  }
  kmer_tail_reference& operator=(kmer_t kmer) {
    kmer_t tail = kmer;
    unsigned char* pos = get_pos();
    for (int i = int(m_ks->m_tail_bytes) - 1; i >= 0; i--) {
      pos[i] = tail & 0xff;
      tail >>= 8;
    }
    return *this;
  }

  bool operator<(const kmer_tail_reference& rhs) const {
    return memcmp(get_pos(), rhs.get_pos(), m_ks->m_tail_bytes) < 0;
  }

  bool operator==(const kmer_tail_reference& rhs) const {
    return memcmp(get_pos(), rhs.get_pos(), m_ks->m_tail_bytes) == 0;
  }

 private:
  unsigned char* get_pos() const {
    return m_ks->m_table + m_table_index * m_ks->m_tail_bytes;
  }
  unsigned get_tail_bytes() const { return m_ks->m_tail_bytes; }
  friend void swap(kmer_tail_reference lhs, kmer_tail_reference rhs) {
    DCHECK_EQ(lhs.m_ks, rhs.m_ks);
    unsigned tail_bytes = lhs.get_tail_bytes();
    char tmp_buf[tail_bytes];
    unsigned char* lhs_pos = lhs.get_pos();
    unsigned char* rhs_pos = rhs.get_pos();
    memcpy(tmp_buf, lhs_pos, tail_bytes);
    memcpy(lhs_pos, rhs_pos, tail_bytes);
    memcpy(rhs_pos, tmp_buf, tail_bytes);
  }

  kmer_set* m_ks = nullptr;
  lookup_t m_table_index = 0;
};

struct kmer_set::kmer_tail_and_flags {
  kmer_t tail;
  unsigned flags;

  kmer_tail_and_flags(const kmer_tail_and_flags_reference& ref);
  bool operator<(const kmer_tail_and_flags& rhs) const { return tail < rhs.tail; }
  bool operator<(const kmer_tail_and_flags_reference& rhs) const;
};

struct kmer_set::kmer_tail_and_flags_reference {
  kmer_tail_reference tail;
  flags_reference flags;

  kmer_tail_and_flags_reference(const kmer_tail_reference& tail_arg,
                                const flags_reference& flags_arg)
      : tail(tail_arg), flags(flags_arg) {}
  kmer_tail_and_flags_reference& operator=(const kmer_tail_and_flags& rhs) {
    tail = rhs.tail;
    flags = rhs.flags;
    return *this;
  }
  kmer_tail_and_flags_reference& operator=(const kmer_tail_and_flags_reference& rhs) {
    tail = rhs.tail;
    flags = rhs.flags;
    return *this;
  }
  friend void swap(kmer_tail_and_flags_reference lhs, kmer_tail_and_flags_reference rhs) {
    swap(lhs.tail, rhs.tail);
    std::swap(lhs.flags, rhs.flags);
  }
  bool operator<(const kmer_tail_and_flags_reference& rhs) const { return tail < rhs.tail; }
  bool operator<(const kmer_tail_and_flags& rhs) const { return tail < rhs.tail; }
};

kmer_set::kmer_tail_and_flags::kmer_tail_and_flags(const kmer_tail_and_flags_reference& ref)
    : tail(ref.tail), flags(ref.flags) {}

bool kmer_set::kmer_tail_and_flags::operator<(const kmer_tail_and_flags_reference& rhs) const {
  return tail < rhs.tail;
}

class kmer_set::kmer_tail_iterator
    : public boost::iterator_facade<kmer_tail_iterator, kmer_tail_and_flags,
                                    std::random_access_iterator_tag,
                                    kmer_tail_and_flags_reference, ptrdiff_t> {
 public:
  kmer_tail_iterator() = default;
  kmer_tail_iterator(kmer_set* ks, lookup_t index)
      : m_ks(ks), m_table_index(index) {}
  kmer_tail_iterator(const kmer_tail_iterator&) = default;
  kmer_tail_iterator& operator=(const kmer_tail_iterator&) = default;

 private:
  friend class boost::iterator_core_access;
  kmer_tail_and_flags_reference dereference() const {
    DCHECK(m_ks);
    return kmer_tail_and_flags_reference(kmer_tail_reference(m_ks, m_table_index),
                                         (*m_ks->m_flags_table)[m_table_index]);
  }
  bool equal(const kmer_tail_iterator& rhs) const {
    return m_table_index == rhs.m_table_index;
  }
  void increment() {
    ++m_table_index;
  }
  void decrement() {
    --m_table_index;
  }
  void advance(ptrdiff_t n) {
    m_table_index += n;
  }
  ptrdiff_t distance_to(const kmer_tail_iterator& rhs) const {
    return ptrdiff_t(rhs.m_table_index) - ptrdiff_t(m_table_index);
  }

  kmer_set* m_ks = nullptr;
  lookup_t m_table_index = std::numeric_limits<lookup_t>::max();
};

constexpr size_t kmer_set::k_not_present;

void kmer_set::kmer_serialized::validate() {
  // Apparently boost's random_access_traversal_tag isn't understood by
  // STL, so make sure STL understands this is a random access iterator
  // so it doesn't spend all its time doing things like incrementing it
  // to measure distance.
  static_assert(std::is_same<std::iterator_traits<kmer_tail_iterator>::iterator_category,
                             std::random_access_iterator_tag>::value,
                "kmer tail iterator must be random access!");

  SPLOG_P(LOG_DEBUG,
          "kmer_set::kmer_serialized::validate> size: %lu, kmer_size: %lu, table: %lu, lookup: %lu",
          size, kmer_size, table.get_size(), lookup.get_size());
}

kmer_t kmer_set::kmer_from_parts(size_t lookup, kmer_t tail) const {
  DCHECK_LE(m_lookup_bits, 2 * m_kmer_size);
  kmer_t shifted_lookup = lookup << (2 * m_kmer_size - m_lookup_bits);
  kmer_t result = shifted_lookup | tail;

  // If there are bits in common, double check they match between lookup and tail.
  unsigned total_bits = m_lookup_bits + 8 * m_tail_bytes;
  DCHECK_GE(total_bits, 2 * m_kmer_size)
      << "Negative overlap between lookup and tail?";
  unsigned overlap_bits = total_bits - 2 * m_kmer_size;
  if (overlap_bits) {
    kmer_t overlap_mask = (kmer_t(1) << overlap_bits) - 1;
    unsigned tail_shift = m_tail_bytes * 8 - overlap_bits;
    DCHECK_EQ(tail_shift, m_kmer_size * 2 - m_lookup_bits);
    DCHECK_EQ(lookup & overlap_mask,
              (tail >> tail_shift) & overlap_mask)
        << overlap_bits << " bits of overlap does not match tail.  lookup=0x" << std::hex << lookup
        << " tail=0x" << tail << std::dec << " tail shift=" << tail_shift;
  }

  return result;
}

kmer_t kmer_set::const_iterator::dereference() const
{
  DCHECK_LT(m_table_index, m_self->m_size);
  DCHECK_LT(m_lookup_index, m_self->m_lookup_size);
  size_t tail_part = m_self->kmer_tail(m_table_index);
  return m_self->kmer_from_parts(m_lookup_index, tail_part);
}

unsigned kmer_set::const_iterator::get_flags() const {
  return m_self->m_flags_table->get(m_table_index);
}

void kmer_set::const_iterator::seek_fixup() {
  auto lookup_pos = std::upper_bound(
      m_self->m_lookup, m_self->m_lookup + m_self->m_lookup_size + 1,
      m_table_index);
  m_lookup_index = lookup_pos - m_self->m_lookup;
  if (m_lookup_index) {
    --m_lookup_index;
  }
}

kmer_set::kmer_set(
	kv_source& source,
	size_t count,
	size_t kmer_size,
	const callback_t& callback)
{
	std::string key;
	std::string value;
	size_t cur = 0;
    kmer_tail_iterator cur_it;
	kmer_t cur_head = 0;
	while (source.read(key, value)) {
		kmer_t all;
		msgpack_deserialize(all, key);
		if (cur == 0) {
			create_sizes(count, kmer_size);
			alloc_tables();
            cur_it = kmer_tail_iterator(this, 0);
		}
		kmer_t head = lookup_for_kmer(all);
		kmer_t tail = tail_for_kmer(all);
		while (cur_head != head) {
			cur_head++;
			m_lookup[cur_head] = lookup_t(cur);
		}
        cur_it->tail = tail;
		callback(cur, all, m_kmer_size, value);
		cur++;
        cur_it++;
	}
	if (cur != count) {
		throw io_exception("Record count incorrect");
	}
	else if (cur == 0) {
		throw io_exception("There are no k-mers in your data.  Try reducing the k-mer filtering minimum score.");
	}
	while (cur_head != m_lookup_size) {
		cur_head++;
		m_lookup[cur_head] = lookup_t(cur);
	}
	m_lookup[m_lookup_size + size_t(1)] = lookup_t(cur + 1);
}

kmer_set::kmer_set(const std::string& serialized, const progress_handler_t& progress)
{
	kmer_serialized ks;
	json_deserialize(ks, serialized);
	ks.validate();
	create_sizes(ks.orig_size, ks.kmer_size);
    m_size = ks.size;
	resource_manager resmgr;
	resmgr.read_resource(m_lookup_buf, ks.lookup, subprogress(progress, 0.0, 0.4));
	resmgr.read_resource(m_table_buf, ks.table, subprogress(progress, 0.4, 1.0));
	m_lookup = (uint32_t *) m_lookup_buf.buffer();
	m_table = (unsigned char*) m_table_buf.buffer();
}

std::string kmer_set::save(const path& root, const progress_handler_t& progress)
{
	kmer_serialized ks;
	resource_manager resmgr;
    ks.orig_size = m_orig_size;
	ks.size = m_size;
	ks.kmer_size = m_kmer_size;
	resmgr.write_resource(ks.lookup, m_lookup_buf, root, "lookup", subprogress(progress, 0.0, 0.4));
	resmgr.write_resource(ks.table, m_table_buf, root, "table", subprogress(progress, 0.4, 1.0));
	ks.validate();
	return json_serialize(ks);
}

kmer_set::const_iterator kmer_set::find(kmer_t key) const {
  size_t table_pos = find_table_index(key);
  if (table_pos == k_not_present) {
    return end();
  }
  return const_iterator(*this, lookup_for_kmer(key), table_pos);
}

template<size_t tail_bytes>
size_t kmer_set::sized_find_internal(kmer_t key) const {
  size_t index = lookup_for_kmer(key);
  kmer_t tail = tail_for_kmer(key);
  // SPLOG("Index = %lu, tail = %s", index, dna_sequence(tail, m_tail_bases).as_string().c_str());
  // SPLOG("Start = %d, end = %u", m_lookup[index], m_lookup[index+1]);

  unsigned char tail_search[tail_bytes];
  for (size_t i = 0; i < tail_bytes; ++i) {
    tail_search[tail_bytes - i - 1] = tail & 0xFF;
    tail >>= 8;
  }

  DCHECK_GT(tail_bytes, 0);
  DCHECK_LT(tail_bytes, sizeof(kmer_t));

  size_t region_start = m_lookup[index];
  size_t region_end = m_lookup[index + 1];
  while (region_start < region_end) {
    size_t pos = (region_start + region_end) / 2;
    unsigned char* ptr = m_table + pos * tail_bytes;

    bool found_it = true;
    for (unsigned i = 0; i != tail_bytes; ++i) {
      if (ptr[i] != tail_search[i]) {
        found_it = false;
        if (ptr[i] < tail_search[i]) {
          region_start = pos + 1;
        } else {  // ptr[i] > tail_search[i]
          region_end = pos;
        }
        break;
      }
    }

    if (found_it) {
      // All bytes are equal.
      return pos;
    }
  }

  // Estimate location of kmer within region
  return k_not_present;
}

size_t kmer_set::find_table_index(kmer_t key) const {
  switch (m_tail_bytes) {
    case 1:
      return sized_find_internal<1>(key);
    case 2:
      return sized_find_internal<2>(key);
    case 3:
      return sized_find_internal<3>(key);
    case 4:
      return sized_find_internal<4>(key);
    case 5:
      return sized_find_internal<5>(key);
    case 6:
      return sized_find_internal<6>(key);
    case 7:
      return sized_find_internal<7>(key);
  }
  LOG(FATAL) << "Invalid number of tail bytes: " << m_tail_bytes;
  return k_not_present;
}

kmer_t kmer_set::kmer_tail(size_t index) const
{
	kmer_t ret = 0;
	for (size_t i = 0; i < m_tail_bytes; i++) {
		ret <<= 8;
		ret |= kmer_t(m_table[m_tail_bytes * index + i]);
	}
	return ret;
}

unsigned kmer_set::get_flags(size_t index) const {
  return m_flags_table->get(index);
}

void kmer_set::create_sizes(size_t size, size_t kmer_size)
{
	m_size = m_orig_size = size;
	m_kmer_size = kmer_size;
	m_tail_bytes = ((kmer_size + 3) / 4); // Start with everything in the tail
	size_t best_memory = m_tail_bytes * size; // Compute initial memory size
	//SPLOG("Computing table for size = %lu, kmer_size = %lu", m_size, m_kmer_size);
	//SPLOG("Also, size_of(size_t) = %lu", sizeof(size_t));
	//SPLOG("Initial guess: tail bytes = %lu, memory = %lu", m_tail_bytes, best_memory);
	while (m_tail_bytes > 0) {
		size_t tail_bytes = m_tail_bytes - 1;
		size_t head_bits = 2 * (kmer_size - (4 * tail_bytes));
		size_t new_memory = (1UL << head_bits) * sizeof(lookup_t) + size * tail_bytes;
		//SPLOG("Next guess: tail bytes = %lu, head_bits = %lu, memory = %lu", tail_bytes, head_bits, new_memory);
		if (new_memory > best_memory) {
			break;
		}
		best_memory = new_memory;
		m_tail_bytes = tail_bytes;
	}
	m_tail_bases = std::min(m_kmer_size, m_tail_bytes * 4);
	m_lookup_bits = 2*(m_kmer_size - m_tail_bases);
	m_lookup_size = (size_t(1) << m_lookup_bits);
}

void kmer_set::create_sizes_with_ram(size_t size, size_t kmer_size, size_t max_ram_bytes) {
  m_size = m_orig_size = size;
  m_kmer_size = kmer_size;
  m_tail_bytes = (kmer_size * 2 + 7) / 8;
  m_lookup_bits = 0;
  m_lookup_size = 1;
  m_tail_bases = std::min(m_kmer_size, m_tail_bytes * 4);

  // Minimum load factor for lookup table to limit sparseness.
  constexpr double k_min_load_factor = 1;
  size_t best_memory = m_tail_bytes * m_size + (k_flag_bits * m_size) / 8;  // Compute initial memory size
  SPLOG("Computing table for size = %lu, kmer_size = %lu", m_size, m_kmer_size);
  SPLOG("Also, size_of(size_t) = %lu", sizeof(size_t));
  SPLOG("Initial guess: tail bytes = %lu, memory = %lu", m_tail_bytes, best_memory);

  // Optimize for lowest number of tail bytes we can use.
  while (m_tail_bytes > 0) {
    size_t tail_bytes = m_tail_bytes - 1;
    size_t lookup_bits = 2 * (kmer_size - (4 * tail_bytes));
    size_t lookup_size = 1UL << lookup_bits;
    size_t new_memory = lookup_size * sizeof(lookup_t) + size * tail_bytes + (k_flag_bits * size) / 8;
    SPLOG("Next guess: tail bytes = %lu, head_bits = %lu, memory = %lu", tail_bytes, lookup_bits,
          new_memory);
    if (new_memory > best_memory) {
      if (new_memory > max_ram_bytes) {
        // Uses too much ram
        SPLOG("Uses more ram than %lu allowed; done search for tail bytes", max_ram_bytes);
        break;
      }

      if (lookup_size * k_min_load_factor > size) {
        // Too sparse
        SPLOG("Load factor %f too sparse; done search for tail bytes", size*1./lookup_size);
        break;
      }

      // Otherwise, continue on to try to get a better one.
    }
    best_memory = new_memory;
    m_tail_bytes = tail_bytes;
    m_lookup_size = lookup_size;
    m_lookup_bits = lookup_bits;
    m_tail_bases = std::min(m_kmer_size, m_tail_bytes * 4);
  }

  // See if we can use more head bits to decrease the amount of binary search we need to do.
  while (m_lookup_bits < kmer_size * 2) {
    size_t new_lookup_bits = m_lookup_bits + 1;
    size_t new_lookup_size = (1UL << new_lookup_bits);
    size_t new_memory = new_lookup_size * sizeof(lookup_t) + size * m_tail_bytes;

    SPLOG("Maybe expand lookup bits: new_lookup_bits = %lu, new_memory = %lu", new_lookup_bits,
          new_memory);

    if (new_memory > max_ram_bytes) {
      SPLOG("Uses more ram than %lu allowed; done search for lookup expansion",
            max_ram_bytes);
      // Uses too much ram
      break;
    }

    if (new_lookup_size * k_min_load_factor > size) {
      // Too sparse
      SPLOG("Load factor %f too sparse; done search for lookup expansion",
            size * 1. / new_lookup_size);
      break;
    }

    m_lookup_size = new_lookup_size;
    m_lookup_bits = new_lookup_bits;
  }

  size_t lookup_mem = m_lookup_size * sizeof(lookup_t);
  size_t tail_mem = size * m_tail_bytes;
  SPLOG(
      "kmer_set: Using %lu bits of lookup and %lu tail bytes; lookup size=%lu MB, tail size=%lu MB "
      "load factor=%f",
      m_lookup_bits, m_tail_bytes, lookup_mem / (1024 * 1024), tail_mem / (1024 * 1024),
      m_size * 1. / m_lookup_size);
  CHECK_GE(m_lookup_bits / 2 + (m_tail_bytes * 4), m_kmer_size);
}

void kmer_set::alloc_tables()
{
	SPLOG("kmer_set> Allocating lookup of %lu, table = %lu", m_lookup_size+size_t(2), m_tail_bytes * m_size);
	SPLOG("kmer_set> m_lookup_bits = %lu, m_tail_bases = %lu", m_lookup_bits, m_tail_bases);
	resource_manager resmgr;
	resmgr.create_resource(m_lookup_buf, sizeof(uint32_t) * (m_lookup_size + 2UL));
	m_lookup = (uint32_t*) m_lookup_buf.buffer();
	resmgr.create_resource(m_table_buf, sizeof(unsigned char) * m_tail_bytes * m_size);
	m_table = (unsigned char*) m_table_buf.buffer();
	m_lookup[0] = 0;
    m_flags_table = make_unique<flags_table_t>(m_size, "kmer_set:flags_table");
}

void kmer_set::alloc_tables_in_memory() {
  SPLOG("kmer_set> Allocating in RAM lookup of %lu, table = %lu",
        m_lookup_size + size_t(2), m_tail_bytes * m_size);
  SPLOG("kmer_set> m_lookup_bits = %lu, m_tail_bases = %lu", m_lookup_bits,
        m_tail_bases);
  m_lookup_membuf = mutable_membuf(new owned_membuf(
      sizeof(uint32_t) * (m_lookup_size + 2UL), "build_kmer_set_lookup"));
  m_lookup = (uint32_t*)m_lookup_membuf.mutable_data();
  m_table_membuf = mutable_membuf(new owned_membuf(
      sizeof(unsigned char) * m_tail_bytes * m_size, "build_kmer_set_table"));
  m_table = (unsigned char*)m_table_membuf.mutable_data();
  m_lookup[0] = 0;
  m_flags_table = make_unique<flags_table_t>(m_size, "kmer_set:flags_table");
}

void kmer_set::save_memory_tables() {
  SPLOG("kmer_set> Saving kmer set to resource manager");
  alloc_tables();
  memcpy(m_lookup_buf.buffer(), m_lookup_membuf.data(),
          sizeof(uint32_t) * (m_lookup_size + 2UL));
  m_lookup_membuf = mutable_membuf();
  memcpy(m_table_buf.buffer(), m_table_membuf.data(),
          sizeof(unsigned char) * m_tail_bytes * m_size);
  m_table_membuf = mutable_membuf();
}

kmer_set::kmer_set(size_t max_count, size_t kmer_size,
                   size_t max_ram,
                   const kmer_source_f& get_kmers,
                   progress_handler_t progress) {
  create_sizes_with_ram(max_count, kmer_size, max_ram);
  alloc_tables_in_memory();
  SPLOG("Generating kmer set for %ld kmers of size %ld", max_count, kmer_size);
  track_mem::reset_stats();

  bool did_limit_size = false;

  if (m_size >= std::numeric_limits<lookup_t>::max()) {
    m_size = std::numeric_limits<lookup_t>::max() - 1;
    SPLOG("Limiting kmer set size build table to %ld", m_size);
    did_limit_size = true;
  }

  // Start out with the lookup table like this:
  // >0 0 0 0 0 0
  // (where ">" specifies where m_lookup points)
  m_lookup++;
  // Count the number of kmers for each prefix into the "lookup" table.
  // Lookup table is now like this:
  // 0 >0 0 0 0 0
  get_kmers(
      [&](kmer_t all, unsigned /* flags */) {
        kmer_t head = lookup_for_kmer(all);
        CHECK_LT(head, m_lookup_size);
        lookup_t* lookup = &m_lookup[head];
        size_t new_pos = __sync_add_and_fetch(lookup, 1);
        if (new_pos > m_size) {
          if (did_limit_size) {
            SPLOG("Too many kmers for kmer table!  Try increasing --min-kmer-count");
            throw(io_exception(
                printstring("Too many kmers for kmer table!  Try increasing --min-kmer-count")));
          } else {
            CHECK_LE(new_pos, m_size)
                << "Overflow of kmer set table size";
          }
        }
      },
      subprogress(progress, 0, 0.4));

  // Lookup table is now like this:
  // 0 >5 8 1 100 2
  // Change from counts to offsets into the tails table.
  size_t tot_tails = 0;
  for (size_t idx = 0; idx != m_lookup_size; ++idx) {
    lookup_t kmer_count = m_lookup[idx];
    CHECK_LT(tot_tails, std::numeric_limits<lookup_t>::max()) << "Too many kmers for kmer table!";
    m_lookup[idx] = tot_tails;
    tot_tails += kmer_count;
  }

  // Lookup table is now like this:
  // 0 >0 5 13 14 114

  // Discard any table entries we didn't need.
  CHECK_LE(tot_tails, m_size);
  SPLOG("After filtering, keeping %lu/%lu kmers (%.2f%%)", tot_tails, m_size,
        tot_tails * 100. / m_size);
  m_size = tot_tails;
  CHECK_GE(m_orig_size, m_size);

  // For each prefix, populate the tails table in no particular order.
  get_kmers(
      [&](kmer_t all, unsigned flags) {
        kmer_t head = lookup_for_kmer(all);
        kmer_t tail = tail_for_kmer(all);
        CHECK_LT(head, m_lookup_size);
        lookup_t* lookup = &m_lookup[head];
        lookup_t cur = __sync_fetch_and_add(lookup, 1);
        CHECK_LT(cur, m_size);
        kmer_tail_iterator cur_it(this, cur);
        CHECK_LT(flags, unsigned(1) << k_flag_bits);
        cur_it->tail = tail;
        cur_it->flags = flags;
      },
      subprogress(progress, 0.4, 0.8));

  // Lookup table is now like this:
  // 0 >5 13 14 114 116
  m_lookup--;

  // Lookup table is now like this, and contains the final offsets
  // into the tails table:
  // >0 5 13 14 114 116

  SPLOG("Sorting kmer set");
  // For each prefix, sort the tails table.
  parallel_for(0, m_lookup_size,
               [&](size_t start, size_t limit) {
                 for (size_t idx = start; idx != limit; ++idx) {
                   lookup_t region_begin = m_lookup[idx];
                   lookup_t region_end = m_lookup[idx + 1];
                   std::sort(kmer_tail_iterator(this, region_begin),
                             kmer_tail_iterator(this, region_end));
                 }
               },
               subprogress(progress, 0.8, 1));
  CHECK_EQ(m_lookup[0], 0);
  CHECK_EQ(m_lookup[m_lookup_size], m_size);
  SPLOG(
      "Done with kmer set generation, lookup size %ld, table size %ld, %ld "
      "tail "
      "bytes, %.2f MB total",
      m_lookup_size, m_size, m_tail_bytes,
      (m_lookup_size * sizeof(lookup_t) + m_size * m_tail_bytes) / 1024. /
          1024.);
}


void kmer_set::copy_into_ram() {
  SPLOG("Loading kmer set into RAM");
  m_lookup_membuf = mutable_membuf(new owned_membuf(
      m_lookup_buf.data(), m_lookup_buf.size(), "kmer_set_lookup"));
  m_table_membuf = mutable_membuf(new owned_membuf(
      m_table_buf.data(), m_table_buf.size(), "kmer_set_table"));
  m_lookup = (lookup_t*) m_lookup_membuf.mutable_data();
  m_table = (unsigned char *) m_table_membuf.mutable_data();
}

#include "modules/io/sparse_multi.h"

const product_version sparse_multi::sparse_multi_version{"1.0.0"};

sparse_multi::sparse_multi(const spiral_file_open_state& state) {
  state.enforce_max_version("sparse_multi", sparse_multi_version);

  m_source_to_mid.reset(new bitcount(state.open_subpart("source_to_mid")));
  m_dest_to_mid.reset(new bitcount(state.open_subpart("dest_to_mid")));
}

sparse_multi::sparse_multi(std::unique_ptr<bitcount> source_to_mid,
                           std::unique_ptr<bitcount> dest_to_mid)
    : m_source_to_mid(std::move(source_to_mid)),
      m_dest_to_mid(std::move(dest_to_mid)) {}

std::pair<uint64_t, uint64_t> sparse_multi::lookup(
    uint64_t source_index) const {
  CHECK_LT(source_index, m_source_to_mid->size());

  if (!m_source_to_mid->get(source_index)) {
    return std::make_pair(0, 0);
  }

  uint64_t mid_index = m_source_to_mid->count(source_index);
  uint64_t dest_index_start = m_dest_to_mid->find_count(mid_index);
  uint64_t dest_index_limit = m_dest_to_mid->find_count(mid_index + 1);

  return std::make_pair(dest_index_start, dest_index_limit);
}

uint64_t sparse_multi::lookup_lower_bound(uint64_t source_index) const {
  CHECK_LT(source_index, m_source_to_mid->size());
  uint64_t mid_index = m_source_to_mid->count(source_index);
  uint64_t dest_index_start = m_dest_to_mid->find_count(mid_index);

  return dest_index_start;
}

std::pair<uint64_t, uint64_t> sparse_multi::lookup_range(
    uint64_t source_index_start, uint64_t source_index_limit) const {
  CHECK_LE(source_index_start, source_index_limit);
  CHECK_LE(source_index_start, m_source_to_mid->size());
  uint64_t mid_index_start = m_source_to_mid->count(source_index_start);
  uint64_t dest_index_start = m_dest_to_mid->find_count(mid_index_start);

  DCHECK_LE(source_index_limit, m_source_to_mid->size());
  uint64_t mid_index_limit = m_source_to_mid->count(source_index_limit);
  uint64_t dest_index_limit = m_dest_to_mid->find_count(mid_index_limit);
  return std::make_pair(dest_index_start, dest_index_limit);
}

uint64_t sparse_multi::reverse_lookup(uint64_t dest_index) const {
  CHECK_LT(dest_index, m_dest_to_mid->size());

  uint64_t mid_index = m_dest_to_mid->count(dest_index);
  if (!m_dest_to_mid->get(dest_index)) {
    mid_index--;
  }
  uint64_t source_index = m_source_to_mid->find_count(mid_index);

  return source_index;
}

uint64_t sparse_multi::lookup_dest_to_mid(uint64_t dest_index) const {
  CHECK_LT(dest_index, m_dest_to_mid->size());
  uint64_t mid_index = m_dest_to_mid->count(dest_index);
  if (!m_dest_to_mid->get(dest_index)) {
    mid_index--;
  }
  return mid_index;
}

uint64_t sparse_multi::lookup_mid_to_source(uint64_t mid_index) const {
  uint64_t source_index = m_source_to_mid->find_count(mid_index);
  return source_index;
}

sparse_multi_builder::sparse_multi_builder(
    const spiral_file_create_state& state, uint64_t n_source_elems,
    uint64_t n_dest_elems) {
  state.set_version("sparse_multi", sparse_multi::sparse_multi_version);

  m_source_to_mid.reset(
      new bitcount(state.create_subpart("source_to_mid"), n_source_elems));
  m_dest_to_mid.reset(
      new bitcount(state.create_subpart("dest_to_mid"), n_dest_elems));
}

uint64_t sparse_multi_builder::add(uint64_t source_index) {
  CHECK(m_source_to_mid);
  CHECK(m_dest_to_mid);

  bool set_dest_bit = false;
  if (m_dest_count == 0) {
    // first item
    set_dest_bit = true;
  }

  if (source_index != m_last_source_seen) {
    CHECK_LT(m_last_source_seen, source_index)
        << "Source indexes must not descend during building";
    set_dest_bit = true;
    m_last_source_seen = source_index;
  }

  m_source_to_mid->set(m_last_source_seen, true);
  if (set_dest_bit) {
    m_dest_to_mid->set(m_dest_count, true);
  }

  return m_dest_count++;
}

std::unique_ptr<sparse_multi> sparse_multi_builder::finalize() {
  CHECK(m_source_to_mid);
  CHECK(m_dest_to_mid);
  CHECK_EQ(m_dest_to_mid->size(), m_dest_count);

  m_source_to_mid->finalize();
  m_dest_to_mid->finalize();

  std::unique_ptr<sparse_multi> result(
      new sparse_multi(std::move(m_source_to_mid), std::move(m_dest_to_mid)));

  m_source_to_mid.reset();
  m_dest_to_mid.reset();

  return result;
}

void sparse_multi_builder::build_from_old_format(const char* gross_ids_buf,
                                                 const char* fine_ids_buf) {
  CHECK_GT(m_source_to_mid->size(), 0);

  const uint32_t* gross_ids = reinterpret_cast<const uint32_t*>(gross_ids_buf);
  const uint16_t* fine_ids = reinterpret_cast<const uint16_t*>(fine_ids_buf);

  size_t n_gross_ids = ((m_source_to_mid->size() - 1) >> 16) + 1;
  size_t cur_dst = 0;

  for (auto gross_index = 0ULL; gross_index < n_gross_ids; gross_index++) {
    auto fine_id_offset = gross_ids[gross_index];
    auto fine_id_end = gross_ids[gross_index + 1];
    for (auto fine_index = fine_id_offset; fine_index < fine_id_end;
         fine_index++) {
      uint64_t seqset_entry_id = (gross_index << 16) + fine_ids[fine_index];
      CHECK_EQ(cur_dst, add(seqset_entry_id));
      cur_dst++;
    }
  }
  CHECK_EQ(cur_dst, m_dest_to_mid->size());
}

sparse_multi::iterator& sparse_multi::iterator::seek_to(uint64_t source_index) {
  m_source_index = source_index;
  CHECK_LE(source_index, m_sm->source_elem_count());
  if (m_source_index == m_sm->source_elem_count()) {
    seek_to_end();
    return *this;
  }
  while (!m_sm->m_source_to_mid->get(m_source_index)) {
    m_source_index++;
    if (m_source_index == m_sm->source_elem_count()) {
      seek_to_end();
      return *this;
    }
  }

  uint64_t mid_index = m_sm->m_source_to_mid->count(m_source_index);
  m_dest_index = m_sm->m_dest_to_mid->find_count(mid_index);
  calculate_next_dest();

  return *this;
}

membuf_cachelist sparse_multi::membufs() const {
  membuf_cachelist results;
  if (m_source_to_mid) {
    results += m_source_to_mid->membufs();
  }
  if (m_dest_to_mid) {
    results += m_dest_to_mid->membufs();
  }
  return results;
}

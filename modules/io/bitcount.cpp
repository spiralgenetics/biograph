#include <cassert>

#include "base/base.h"
#include "modules/io/bitcount.h"
#include "modules/io/log.h"
#include "modules/test/coverage.h"

DECLARE_TEST_COVERAGE(bitcount);

const product_version bitcount::bitcount_version{"1.0.0"};

static size_t round_div(size_t size, size_t div) {
  return (size + (div - 1)) / div;
}

size_t bitcount::compute_size(size_t nbits) {
  return bits_mem_size(nbits) + subaccum_mem_size(nbits) +
         accum_mem_size(nbits);
}

bitcount::bitcount(const void *const_buf, size_t nbits)
    : m_nbits(nbits), m_mutable(true) {
  void *buf = const_cast<void *>(const_buf);
  mutable_membuf borrowed(
      new borrowed_mutable_membuf((char *)buf, compute_size(nbits)));
  m_bits = m_mutable_bits = borrowed.subbuf(0, bits_mem_size(nbits));
  m_subaccum = m_mutable_subaccum =
      borrowed.subbuf(bits_mem_size(nbits), subaccum_mem_size(nbits));
  m_accum = m_mutable_accum = borrowed.subbuf(
      bits_mem_size(nbits) + subaccum_mem_size(nbits), accum_mem_size(nbits));
}

bitcount::bitcount(size_t nbits) : m_nbits(nbits), m_mutable(true) {
  m_bits = m_mutable_bits = new owned_membuf(bits_mem_size(nbits), "bitcount");
  m_subaccum = m_mutable_subaccum =
      new owned_membuf(subaccum_mem_size(nbits), "bitcount_subaccum");
  m_accum = m_mutable_accum =
      new owned_membuf(accum_mem_size(nbits), "bitcount_accum");
}

bitcount::bitcount(const spiral_file_create_state &state, size_t nbits)
    : m_nbits(nbits), m_mutable(true) {
  state.set_version("bitcount", bitcount_version);
  bc_metadata md;
  md.nbits = nbits;
  state.create_json<bc_metadata>("bitcount.json", md);
  m_bits = m_mutable_bits = state.create_membuf("bits", bits_mem_size(nbits));
  m_subaccum = m_mutable_subaccum =
      state.create_membuf("subaccum", subaccum_mem_size(nbits));
  m_accum = m_mutable_accum =
      state.create_membuf("accum", accum_mem_size(nbits));
}

bitcount::bitcount(const spiral_file_open_state &state) : m_mutable(false) {
  state.enforce_max_version("bitcount", bitcount_version);
  bc_metadata bc = state.open_json<bc_metadata>("bitcount.json");
  m_nbits = bc.nbits;
  m_bits = state.open_membuf("bits");
  CHECK_EQ(bits_mem_size(m_nbits), m_bits.size());
  m_subaccum = state.open_membuf("subaccum");
  CHECK_EQ(subaccum_mem_size(m_nbits), m_subaccum.size());
  m_accum = state.open_membuf("accum");
  CHECK_EQ(accum_mem_size(m_nbits), m_accum.size());
}

size_t bitcount::bits_mem_size(size_t nbits) {
  return round_div(nbits, 64) * sizeof(uint64_t);
}

size_t bitcount::subaccum_mem_size(size_t nbits) {
  return round_div(nbits, 512) * sizeof(uint64_t);
}

size_t bitcount::accum_mem_size(size_t nbits) {
  // Storage for accum, needs room 1 extra bit
  return round_div(nbits + 1, 512) * sizeof(uint64_t);
}

void bitcount::init() {
  CHECK(m_mutable);
  memset(m_mutable_bits.mutable_data(), 0, bits_mem_size(m_nbits));
}

size_t bitcount::finalize(progress_handler_t prog) {
  CHECK(m_mutable);

  if (size() == 0) {
    mutable_accum()[0] = 0;
    return 0;
  }

  uint64_t subaccum = 0;
  uint64_t total = 0;
  for (size_t i = 0; i < round_div(m_nbits, 64); i++) {
    prog(0.0);
    if (i % 8 == 0) {
      mutable_accum()[i / 8] = total;
      if (i != 0) mutable_subaccum()[i / 8 - 1] = subaccum;
      subaccum = 0;
    }
    subaccum <<= 8;
    uint64_t bits = this->bits()[i];
    uint64_t subtotal =
        __builtin_popcount(bits >> 32) + __builtin_popcount(bits & 0xffffffff);
    subaccum |= subtotal;
    total += subtotal;
  }
  size_t left_over = round_div(m_nbits, 64) % 8;
  while (left_over != 0) {
    subaccum <<= 8;
    left_over++;
    left_over %= 8;
  }
  mutable_subaccum()[round_div(m_nbits, 512) - 1] = subaccum;

  // Handle special case of exact even size
  if (m_nbits % 512 == 0) {
    NOTE_TEST_COVERAGE(bitcount);
    mutable_accum()[m_nbits / 512] = total;
  }

  return total;
}

namespace {

// Returns the bit index of the count-th set bit in val.
unsigned select_bit(uint64_t val, unsigned count) {
  unsigned index = 0;

  for (unsigned i = 32; i > 0; i /= 2) {
    uint64_t mask = (1UL << i) - 1;
    unsigned p = __builtin_popcount(val & mask);
    if (p <= count) {
      index += i;
      val >>= i;
      count -= p;
    }
  }
  return index;
}

}  // namespace

void bitcount::make_find_count_index() {
  m_count_set_to_index.clear();
  size_t max_shifted = total_bits() >> k_find_count_count_bits;
  m_count_set_to_index.resize(max_shifted + 2);
  m_count_set_to_index[0] = 0;
  CHECK_LT(size() >> k_find_count_index_bits,
           std::numeric_limits<uint32_t>::max());
  m_count_set_to_index[max_shifted + 1] = size() >> k_find_count_index_bits;

  size_t cur_shifted = 0;
  size_t cur_count = 0;

  size_t idx64_max = round_div(size(), 64);
  for (size_t idx64 = 0; idx64 < idx64_max; ++idx64) {
    uint64_t cur_val64 = bits()[idx64];
    size_t new_count = cur_count + __builtin_popcountl(cur_val64);
    size_t new_shifted = new_count >> k_find_count_count_bits;
    if (new_shifted == cur_shifted) {
      // Fast path; we don't need to generate any index entries for
      // this "accum" so we can skip right to the next one.
      NOTE_TEST_COVERAGE(bitcount);
      cur_count = new_count;
      continue;
    }

    size_t idx = idx64 * 64;
    size_t end_idx = std::min<uint64_t>(size(), (idx64 + 1) * 64);
    while (idx != end_idx) {
      if (cur_val64 & 1) {
        cur_count++;
        if ((cur_count >> k_find_count_count_bits) != cur_shifted) {
          cur_shifted++;
          CHECK_EQ(cur_shifted, cur_count >> k_find_count_count_bits);
          m_count_set_to_index[cur_shifted] = idx >> k_find_count_index_bits;

          NOTE_TEST_COVERAGE(bitcount);
        }
      }
      cur_val64 >>= 1;
      idx++;
    }
    CHECK_EQ(cur_shifted, new_shifted);
    if (idx < size()) {
      NOTE_TEST_COVERAGE(bitcount);
      CHECK_EQ(new_count, cur_count);
      CHECK_EQ(0, cur_val64);
    } else {
      NOTE_TEST_COVERAGE(bitcount);
    }
  }
  CHECK_EQ(max_shifted, cur_shifted);
  CHECK_EQ(total_bits(), cur_count);
}

size_t bitcount::find_count(size_t target_count) const {
  if (size() == 0) {
    CHECK_EQ(0, target_count);
    return 0;
  }

  size_t accum_start = 0;
  size_t shifted = (target_count + 1) >> k_find_count_count_bits;
  size_t accum_max = round_div(size(), 512);

  if (m_count_set_to_index.empty()) {
    accum_start = std::lower_bound(accum(), accum() + round_div(size(), 512),
                                   target_count + 1) -
                  accum();
    CHECK_GT(accum_start, 0);
    accum_start--;

    NOTE_TEST_COVERAGE_IF(bitcount, (size() % 512) != 0);
    NOTE_TEST_COVERAGE_IF(bitcount, (size() % 512) == 0);
    NOTE_TEST_COVERAGE_IF(bitcount, accum()[accum_start] == 0);
    NOTE_TEST_COVERAGE_IF(bitcount,
                          accum()[accum_start] == round_div(size(), 512));
    NOTE_TEST_COVERAGE_IF(bitcount, accum()[accum_start] == 0);
  } else {
    // Use index to find a much smaller range of accum to search.
    accum_start = m_count_set_to_index[shifted];
    // accum_end = m_count_set_to_index[shifted + 1] + 1;
    NOTE_TEST_COVERAGE_IF(bitcount, (size() & 511) != 0);
    NOTE_TEST_COVERAGE_IF(bitcount, (size() & 511) == 0);
  }

  size_t accum_lb = accum_start;
  int linear_search_counter = 0;
  while (accum_lb < accum_max && accum()[accum_lb] < target_count + 1) {
    linear_search_counter++;
    if (linear_search_counter == k_find_count_max_linear_search) {
      NOTE_TEST_COVERAGE(bitcount);
      // Give up on the linear search, and fall back to binary search
      // which is much less friendly with respect to memory latency.
      CHECK(!m_count_set_to_index.empty())
          << "Without a find count index, we should've already done a binary "
             "search.";
      size_t accum_end = m_count_set_to_index[shifted + 1] + 1;
      if (accum_end > accum_max) {
        accum_end = accum_max;
        NOTE_TEST_COVERAGE(bitcount);
      }
      CHECK_LT(accum_lb, accum_end);
      accum_lb = std::lower_bound(accum() + accum_lb, accum() + accum_end,
                                  target_count + 1) -
                 accum();
      NOTE_TEST_COVERAGE_IF(bitcount, accum_lb == accum_max);
      break;
    }
    accum_lb++;
  }
  NOTE_TEST_COVERAGE_IF(bitcount, accum_lb == accum_max);
  accum_lb--;

  size_t cur_count = accum()[accum_lb];
  size_t subaccum = this->subaccum()[accum_lb];

  for (int i = 0; i < 8; i++) {
    size_t this_subaccum = (subaccum >> (56 - i * 8)) & 0xFF;
    NOTE_TEST_COVERAGE_IF(bitcount, i == 7 && this_subaccum == 1);
    NOTE_TEST_COVERAGE_IF(bitcount, i == 7 && target_count == 1);
    if (this_subaccum + cur_count <= target_count) {
      cur_count += this_subaccum;
      continue;
    }

    size_t idx = (accum_lb * 8 + i) * 64;
    uint64_t dat = bits()[accum_lb * 8 + i];
    NOTE_TEST_COVERAGE_IF(bitcount, dat == 0xFFFFFFFFFFFFFFFFUL);
    return idx + select_bit(dat, target_count - cur_count);
  }

  if (accum_lb == (accum_max - 1)) {
	  NOTE_TEST_COVERAGE(bitcount);
	  return size();
  }
  LOG(FATAL) << "Bit for count " << target_count
             << " was not in expected region; size = " << size()
             << " total bits = " << total_bits();
}

membuf_cachelist bitcount::membufs() const {
  return {m_bits, m_accum, m_subaccum};
}


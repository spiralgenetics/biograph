#include "base/base.h"
#include "modules/bio_base/kmer.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/io/utils.h"

#include <cstring>
#include <cstdio>
#include <byteswap.h>

#include <boost/format.hpp>

dna_sequence::dna_sequence(const dna_sequence& rhs)
	: m_data(rhs.m_size ? new unsigned char[rhs.m_size / 4 + 1] : NULL)
	, m_capacity(rhs.m_size)
	, m_size(rhs.m_size)
{
	if (m_size)
		memcpy(m_data, rhs.m_data, rhs.isize());
}

dna_sequence::dna_sequence(dna_sequence&& rhs)
	: m_data(rhs.m_data)
	, m_capacity(rhs.m_size)
	, m_size(rhs.m_size)
{
	rhs.m_data = nullptr;
	rhs.m_capacity = rhs.m_size = 0;
}

dna_sequence::dna_sequence(const const_iterator& start, const const_iterator& end) :
    dna_sequence(dna_slice(start, end)) {}

dna_sequence::dna_sequence(const dna_slice& rhs) {
  if (rhs.empty()) {
    return;
  }
  resize(rhs.size());
  CHECK(copy_bases(rhs, begin()) == end());
}

dna_sequence::dna_sequence(size_t size)
	: m_data(size ? new unsigned char[size/4 + 1] : NULL)
	, m_capacity(size)
	, m_size(size)
{
	if (size == 0)
		return;
	memset(m_data, 0, isize());
	m_data[0] |= (m_size & 3) << 6;
}

dna_sequence::dna_sequence(const dna_base& b, size_t size)
    : m_data(size ? new unsigned char[size / 4 + 1] : NULL), m_capacity(size), m_size(size) {
  if (!size) {
    return;
  }
  CHECK_GT(isize(), 0);
  m_data[isize() - 1] = 0;
  for (size_t i = 0; i < size; i++) {
    (*this)[i] = b;
  }
  m_data[0] &= 0x3f;
  m_data[0] |= (m_size & 3) << 6;
}

dna_sequence::dna_sequence(kmer_t kmer, int size)
        : m_data(new unsigned char[size / 4 + 1])
        , m_capacity(size)
        , m_size(size)
{
	m_data[isize() - 1] = 0;
	for(int i = size - 1; i >= 0; i--)
	{
		(*this)[i] = dna_base((int) (kmer & 0x3));
		kmer >>= 2;
	}
	m_data[0] &= 0x3f;
	m_data[0] |= (m_size & 3) << 6;
}

dna_sequence& dna_sequence::operator=(const dna_sequence& rhs)
{
	if (m_data != rhs.m_data) {
		delete[] m_data;
		m_size = rhs.m_size;
		m_data = rhs.m_size ? new unsigned char[rhs.isize()] : nullptr;
		m_capacity = m_size;
		if (m_size) {
			memcpy(m_data, rhs.m_data, rhs.isize());
		}
	}
	
	return *this;
}

void dna_sequence::resize(size_t new_size) {
  // If we do support resizing smaller, it will need to be cleared out.
  CHECK_GE(new_size, m_size) << "Resizing a dna sequence smaller is not supported.";
  reserve(new_size);
  m_size = new_size;
  if (m_data) {
	m_data[0] &= 0x3f;
	m_data[0] |= (m_size & 3) << 6;
  }
}

dna_sequence& dna_sequence::operator=(const dna_slice& rhs) {
  clear();
  if (rhs.empty()) {
    return *this;
  }
  resize(rhs.size());
  CHECK(copy_bases(rhs, begin()) == end());
  return *this;
}

dna_sequence& dna_sequence::operator=(dna_sequence&& rhs)
{
	if (m_data != rhs.m_data) {
		std::swap(m_data, rhs.m_data);
		std::swap(m_size, rhs.m_size);
		std::swap(m_capacity, rhs.m_capacity);
	}
	
	return *this;
}

dna_sequence::~dna_sequence()
{
	delete[] m_data;
}

void dna_sequence::clear()
{
	m_size = 0;
}

void dna_sequence::push_back(dna_base b)
{
	if (m_size == m_capacity)
	{
		if (m_size < 64)
		{
			reserve(127);
		}
		else
		{
			reserve(m_capacity * 2);
		}
	}
    *end() = b;
    resize(m_size + 1);
}

void dna_sequence::reserve(unsigned long new_capacity)
{
	if (m_capacity >= new_capacity) return;
	if (new_capacity < 127) m_capacity = 127;
	else m_capacity = new_capacity;

	unsigned char* new_data = new unsigned char[m_capacity / 4 + 1];
	memset(new_data, 0, m_capacity / 4 + 1);
	if (m_size != 0)
		memcpy(new_data, m_data, isize());
	delete[] m_data;
	m_data = new_data;
}

dna_sequence::dna_sequence(const std::string& seq, bool packed)
{
	if (packed)
	{
		if (seq.size() == 0)
		{
			m_data = NULL;
			m_size = 0;
			return;
		}
		m_data = new unsigned char[seq.size()];
		memcpy(m_data, seq.data(), seq.size());
		m_size = (seq.size() - 1) * 4 + (m_data[0] >> 6);
		m_capacity = m_size;
	}
	else
	{
		m_size = seq.size();
		if (m_size == 0) { m_data = NULL; return; }
		m_data = new unsigned char[isize()];
		m_data[m_size/4] = 0;
		for(size_t i = 0; i < size(); i++)
			(*this)[i] = dna_base(seq[i]);
		m_data[0] &= 0x3f;
		m_data[0] |= (m_size & 3) << 6;
		m_capacity = m_size;
	}
}

std::string dna_sequence::as_string() const
{
	std::string out;
	for(size_t i = 0; i < size(); i++)
		out.push_back((char) (*this)[i]);
	return out;
}

std::string dna_slice::as_string() const
{
	std::string out;
	for(size_t i = 0; i < size(); i++)
		out.push_back((char) (*this)[i]);
	return out;
}

std::string dna_sequence::as_packed() const
{
	if (m_size == 0) return std::string("\0", 1);
	CHECK((m_data[0] & 0xc0) == (m_size & 3) << 6);
	return std::string((char*) m_data, isize());
}

kmer_t dna_sequence::as_kmer() const
{
	if (m_size > 32)
		throw io_exception("Maximum k-mer size is 32");

	kmer_t r = 0;
	for(size_t i = 0; i < m_size; i++)
	{
		r <<= 2;
		r |= (int) (*this)[i];
	}
	return r;
}

kmer_t dna_slice::as_kmer() const
{
	if (m_size > 32)
		throw io_exception("Maximum k-mer size is 32");

	kmer_t r = 0;
	for(size_t i = 0; i < m_size; i++)
	{
		r <<= 2;
		r |= (int) (*this)[i];
	}
	return r;
}

bool dna_sequence::operator<(const dna_sequence& rhs) const
{
	return subseq_lessthan(begin(), rhs.begin(), size(), rhs.size());
}

bool dna_sequence::operator==(const dna_sequence& rhs) const
{
	if (rhs.size() != m_size) {
		return false;
	}
	if (m_size == 0) {
		return true;
	}
	return memcmp(m_data, rhs.m_data, isize()) == 0;
}

bool dna_sequence::operator==(const dna_slice& rhs) const
{
  return dna_slice(*this) == rhs;
}

dna_sequence dna_sequence::subseq(size_t offset, size_t len) const
{
	dna_sequence r(len);
	for(size_t i = 0; i < len; i++)
		r[i] = (*this)[offset + i];
	return r;
}

dna_sequence dna_sequence::reverse() const
{
	dna_sequence r(m_size);
	for(size_t i = 0; i < size(); i++)
		r[i] = (*this)[size() - 1 - i];
	return r;
}

dna_sequence dna_sequence::rev_comp() const
{
	dna_sequence r(m_size);
	for(size_t i = 0; i < size(); i++)
		r[i] = (*this)[size() - 1 - i].complement();
	return r;
}

dna_sequence dna_sequence::canonicalize(bool& flipped) const
{
	dna_sequence cs = this->rev_comp();
	if (cs < *this)
	{
		flipped = true;
		return cs;
	}
	else
	{
		flipped = false;
		return *this;
	}
}

dna_sequence& dna_sequence::operator+=(const dna_sequence& rhs)
{
  return (*this) += dna_slice(rhs);
}

dna_sequence& dna_sequence::operator+=(const dna_base& rhs)
{
	push_back(rhs);
	
	return *this;
}

dna_sequence& dna_sequence::operator+=(const dna_slice& rhs)
{
  size_t orig_size = size();
  resize(orig_size + rhs.size());
  CHECK(copy_bases(rhs, begin() + orig_size) == end());
  return *this;
}


kmer_t rev_comp(kmer_t in, uint32_t size)
{
	return long_rev_comp_bases(in) >> (8 * sizeof(kmer_t) - 2 *size);
}

kmer_t canonicalize(kmer_t in, int size)
{
	kmer_t cin = rev_comp(in, size);
	return std::min(in, cin);
}

kmer_t canonicalize(kmer_t in, int size, bool& flipped)
{
	kmer_t cin = rev_comp(in, size);
	if (cin < in)
	{
		flipped = true;
		return cin;
	}
	flipped = false;
	return in;
}

inline uint64_t offset_mask(uint8_t offset)
{
	static_assert((0xffffffffffffffffUL >> 2) == 0x3fffffffffffffffUL, "Unsigned ints right shift is arithmetic (sign extension)!");
	CHECK_LT(offset, 64);
	return 0xffffffffffffffffUL >> offset * 2;
}

inline uint64_t complement_offset_mask(uint8_t offset)
{
	static_assert((0xffffffffffffffffUL >> 2) == 0x3fffffffffffffffUL, "Unsigned ints right shift is arithmetic (sign extension)!");
	CHECK_LT(offset, 64);
	return ~(0xffffffffffffffffUL >> 2 * (offset + 1));
}

inline uint64_t length_mask(uint8_t length)
{
	return (length != 32) ? ~offset_mask(length % 32) : 0xffffffffffffffffUL;
}

uint8_t byte_rev_comp_bases(uint8_t a_byte)
{
	uint8_t temp = ((a_byte >> 2) & 0x33) | ((a_byte & 0x33) << 2);
	return ~(((temp >> 4) & 0x0f) | ((temp & 0x0f) << 4));
}

uint64_t long_rev_comp_bases(uint64_t a_byte)
{
	uint64_t temp = ((a_byte >> 2) & 0x3333333333333333UL) | ((a_byte & 0x3333333333333333UL) << 2);
	temp = ((temp >> 4) & 0x0f0f0f0f0f0f0f0fUL) | ((temp & 0x0f0f0f0f0f0f0f0fUL) << 4);
	temp = ((temp >> 8) & 0x00ff00ff00ff00ffUL) | ((temp & 0x00ff00ff00ff00ffUL) << 8);
	temp = ((temp >> 16) & 0x0000ffff0000ffffUL) | ((temp & 0x0000ffff0000ffffUL) << 16);
	return ~((temp >> 32) | (temp << 32));
}

std::ostream& operator<<(std::ostream& os, dna_compare_result compare_result) {
  switch(compare_result) {
    case dna_compare_result::FIRST_IS_LESS: return os << "FIRST_IS_LESS";
    case dna_compare_result::FIRST_IS_PREFIX: return os << "FIRST_IS_PREFIX";
    case dna_compare_result::EQUAL: return os << "EQUAL";
    case dna_compare_result::SECOND_IS_PREFIX: return os << "SECOND_IS_PREFIX";
    case dna_compare_result::SECOND_IS_LESS: return os << "SECOND_IS_LESS";
    default:
      return os << "Unknown compare result " << int(compare_result);
  }
}


namespace {

uint64_t get_fwd_compare_block(const uint8_t* ptr, unsigned offset,
                               unsigned copy_bases) {
  uint64_t result = 0;
  unsigned copy_bytes = (copy_bases + offset + 3) / 4;

  DCHECK_LE(copy_bytes, 8);
  DCHECK_GT(copy_bytes, 0);
  if (copy_bytes > 8 || copy_bytes == 0) {
    // This lets some versions of GCC inline the memcpy call below.
    __builtin_unreachable();
  }
  memcpy(&result, ptr, copy_bytes);
  return result;
}

uint64_t get_rc_compare_block(const uint8_t* ptr, unsigned offset,
                              unsigned copy_bases) {
  uint64_t result = 0;
  unsigned copy_bytes = (copy_bases + 3 + 3 - offset) / 4;
  uint8_t* result_ptr = reinterpret_cast<uint8_t*>(&result);

  DCHECK_LE(copy_bytes, 8);
  DCHECK_GT(copy_bytes, 0);
  if (copy_bytes > 8 || copy_bytes == 0) {
    // This lets some versions of GCC inline the memcpy call below.
    __builtin_unreachable();
  }
  memcpy(result_ptr + sizeof(uint64_t) - copy_bytes, ptr + 1 - copy_bytes,
         copy_bytes);
  return result;
}

// Declare this as ATTRIBUTE_NO_SANITIZE_ADDRESS since the last byte
// sometimes is off the end of the sequence.
uint64_t ATTRIBUTE_NO_SANITIZE_ADDRESS get_full_block(const uint8_t* ptr,
                                                      unsigned offset,
                                                      bool rc) {
  if (rc) {
    return *reinterpret_cast<const uint64_t*>(ptr + 1 - sizeof(uint64_t));
  } else {
    return *reinterpret_cast<const uint64_t*>(ptr);
  }
}

// be64toh_and_rc(val) == long_rev_comp_bases(be64toh(val)), but does a lot
// less work.
uint64_t be64toh_and_rc(uint64_t val) {
  DCHECK_EQ(be64toh(0x0123456789abcdefULL), 0xefcdab8967452301ULL)
      << "This routine will need to be reworked to function on big-endian "
         "architectures.";

  val =
      ((val >> 2) & 0x3333333333333333UL) | ((val & 0x3333333333333333UL) << 2);
  val =
      ((val >> 4) & 0x0f0f0f0f0f0f0f0fUL) | ((val & 0x0f0f0f0f0f0f0f0fUL) << 4);
  return ~val;
}

template <bool compare_full_block>
inline int64_t compare_shifted(const uint8_t*& lhs, unsigned lhs_offset,
                               bool lhs_rc, const uint8_t*& rhs,
                               unsigned rhs_offset, bool rhs_rc,
                               unsigned compare_size) {
  DCHECK_LE(lhs_offset, 3);
  DCHECK_LE(rhs_offset, 3);
  DCHECK_LE(compare_size, 28);
  if (compare_full_block) {
    DCHECK_EQ(compare_size, 28);
  }

  uint64_t lhs_block, rhs_block;
  if (compare_full_block) {
    lhs_block = get_full_block(lhs, lhs_offset, lhs_rc);
    rhs_block = get_full_block(rhs, rhs_offset, rhs_rc);
  } else {
    if (lhs_rc) {
      lhs_block = get_rc_compare_block(lhs, lhs_offset, compare_size);
    } else {
      lhs_block = get_fwd_compare_block(lhs, lhs_offset, compare_size);
    }
    if (rhs_rc) {
      rhs_block = get_rc_compare_block(rhs, rhs_offset, compare_size);
    } else {
      rhs_block = get_fwd_compare_block(rhs, rhs_offset, compare_size);
    }
  }

  if (lhs_rc) {
    lhs_block = be64toh_and_rc(lhs_block);
  } else {
    lhs_block = be64toh(lhs_block);
  }

  if (rhs_rc) {
    rhs_block = be64toh_and_rc(rhs_block);
  } else {
    rhs_block = be64toh(rhs_block);
  }

  uint64_t mask = (1ULL << (64 - 8)) - 1;
  if (!compare_full_block) {
    mask &= (~0ULL) << 2 * (7 * 4 - compare_size);
  }

  lhs_block >>= 8 - 2 * (lhs_rc ? (3 - lhs_offset) : lhs_offset);
  rhs_block >>= 8 - 2 * (rhs_rc ? (3 - rhs_offset) : rhs_offset);

  int64_t res = int64_t(lhs_block & mask) - int64_t(rhs_block & mask);
  if (res) {
    return res;
  }

  if (lhs_rc) {
    lhs -= 7;
  } else {
    lhs += 7;
  }
  if (rhs_rc) {
    rhs -= 7;
  } else {
    rhs += 7;
  }

  return 0;
}

dna_compare_result compare_internal(const uint8_t* lhs, unsigned lhs_offset,
                                    bool lhs_rc, size_t lhs_size,
                                    const uint8_t* rhs, unsigned rhs_offset,
                                    bool rhs_rc, size_t rhs_size) {
  size_t compare_size_left = std::min(lhs_size, rhs_size);

  while (compare_size_left > 7 * 4) {
    auto res = compare_shifted<true /* full block */>(
        lhs, lhs_offset, lhs_rc, rhs, rhs_offset, rhs_rc, 7 * 4);
    if (res) {
      if (res < 0) {
        return dna_compare_result::FIRST_IS_LESS;
      } else {
        return dna_compare_result::SECOND_IS_LESS;
      }
    }
    compare_size_left -= 7 * 4;
  }

  if (compare_size_left) {
    auto res = compare_shifted<false /* not a full block */>(
        lhs, lhs_offset, lhs_rc, rhs, rhs_offset, rhs_rc, compare_size_left);
    if (res) {
      if (res < 0) {
        return dna_compare_result::FIRST_IS_LESS;
      } else {
        return dna_compare_result::SECOND_IS_LESS;
      }
    }
  }

  if (lhs_size == rhs_size) {
    return dna_compare_result::EQUAL;
  } else if (lhs_size < rhs_size) {
    return dna_compare_result::FIRST_IS_PREFIX;
  } else {
    return dna_compare_result::SECOND_IS_PREFIX;
  }
}

}  // namespace

dna_compare_result type_decode<true>::do_compare(data_type data1, ptrdiff_t offset1, bool rev_comp1,
                                                 data_type data2, ptrdiff_t offset2, bool rev_comp2,
                                                 size_t len1, size_t len2) {
  if (len1 == 0) {
    if (len2 == 0) {
      return dna_compare_result::EQUAL;
    } else {
      return dna_compare_result::FIRST_IS_PREFIX;
    }
  } else if (len2 == 0) {
    return dna_compare_result::SECOND_IS_PREFIX;
  }

  CHECK(data1);
  CHECK(data2);
  CHECK_LT(offset1, 4);
  CHECK_LT(offset2, 4);

  const uint8_t* lhs = (const uint8_t*)data1;
  const uint8_t* rhs = (const uint8_t*)data2;

  return compare_internal(lhs, offset1, rev_comp1, len1, rhs, offset2, rev_comp2, len2);
}

unsigned dna_slice::shared_prefix_length(const dna_slice& rhs_slice) const {
  auto lhs = m_begin.get_data_byte();
  auto rhs = rhs_slice.m_begin.get_data_byte();
  auto lhs_offset = m_begin.get_offset_in_byte();
  auto rhs_offset = rhs_slice.m_begin.get_offset_in_byte();
  bool lhs_rc = m_begin.is_rev_comp();
  bool rhs_rc = rhs_slice.m_begin.is_rev_comp();
  auto lhs_size = size();
  auto rhs_size = rhs_slice.size();

  if (!lhs_size || !rhs_size) {
    return 0;
  }

  CHECK_LT(lhs_offset, 4);
  CHECK_LT(rhs_offset, 4);

  unsigned same_bases = 0;

  size_t compare_size_left = std::min(lhs_size, rhs_size);
  while (compare_size_left >= 7 * 4) {
    auto res = compare_shifted<true /* full block */>(
        lhs, lhs_offset, lhs_rc, rhs, rhs_offset, rhs_rc, 7 * 4);
    if (res != 0) {
      break;
    }
    compare_size_left -= 7 * 4;
    same_bases += 7 * 4;
  }

  dna_const_iterator lhs_it = m_begin + same_bases;
  dna_const_iterator rhs_it = rhs_slice.m_begin + same_bases;

  while (compare_size_left && *lhs_it == *rhs_it) {
    ++same_bases;
    --compare_size_left;
    ++lhs_it;
    ++rhs_it;
  }

  return same_bases;
}

dna_sequence::iterator dna_sequence::copy_bases(const dna_slice& slice, iterator dest) {
  CHECK(!dest.is_rev_comp())
      << "copy_bases may only be used with a forward facing iterator as a destination";
  auto it = slice.begin();
  auto end_it = slice.end();

  if (it == end_it) {
    return dest;
  }

  // Advance until we're on a byte boundary in the destination.
  while (dest.get_offset_in_byte() != 0) {
    *dest = *it;
    ++it;
    ++dest;
    if (it == end_it) {
      return dest;
    }
  }

  const uint8_t* src_ptr = it.get_data_byte();
  unsigned src_offset = it.get_offset_in_byte();
  bool src_rc = it.is_rev_comp();
  uint8_t* dest_ptr = dest.get_data_byte();
  unsigned copy_bytes = (end_it - it) / 4;

  if (src_offset == 0 && !src_rc) {
    // No need shift or reverse!  Memcpy will suffice for this case.
    memcpy(dest_ptr, src_ptr, copy_bytes);
    dest_ptr += copy_bytes;
    src_ptr += copy_bytes;
    dest += copy_bytes * 4;
    it += copy_bytes * 4;
  } else {
    if (copy_bytes > 0) {
      // Make sure we don't go off the end of any buffers by a byte,
      // since we're copying in blocks of 8 bytes instead of 7.
      --copy_bytes;
    }
    unsigned fast_iters = copy_bytes / 7;
    unsigned block_shift = 2 * (src_rc ? (3 - src_offset) : src_offset);
    for (unsigned i = 0; i < fast_iters; ++i) {
      uint64_t block = get_full_block(src_ptr, src_offset, src_rc);

      if (src_rc) {
        block = be64toh_and_rc(block);
      } else {
        block = be64toh(block);
      }
      block <<= block_shift;
      *reinterpret_cast<uint64_t*>(dest_ptr) = htobe64(block);

      dest_ptr += 7;
      if (src_rc) {
        src_ptr -= 7;
      } else {
        src_ptr += 7;
      }
    }

    dest += fast_iters * 7 * 4;
    it += fast_iters * 7 * 4;
  }
  CHECK_EQ((void*)src_ptr, (void*)it.get_data_byte());
  CHECK_EQ((void*)dest_ptr, (void*)dest.get_data_byte());
  if (it != end_it) {
    CHECK_EQ(0, dest.get_offset_in_byte());
  }

  while (it != end_it) {
    *dest = *it;
    ++it;
    ++dest;
  }
  return dest;
}

namespace dna_testutil {
std::ostream& default_dna_printer(std::ostream& os, const dna_slice& seq) {
  return os << seq.as_string();
}

std::function<std::ostream&(std::ostream&, const dna_slice&)> g_dna_printer =
    default_dna_printer;
}  // namespace dna_testutil

std::ostream& operator<<(std::ostream& os, const dna_sequence& seq) {
  dna_slice slice(seq.begin(), seq.size());
  return dna_testutil::g_dna_printer(os, slice);
}

std::ostream& operator<<(std::ostream&os, const dna_slice&slice) {
  return dna_testutil::g_dna_printer(os, slice);
}

#pragma once

#include "modules/bio_base/dna_sequence.h"
#include <iterator>

struct kcount_pair
{
	kcount_pair() :
		fwd(0), rev(0) 
	{}
	kcount_pair(uint32_t _fwd, uint32_t _rev) :
		fwd(_fwd), rev(_rev) 
	{}
	
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(fwd);
		FIELD(rev);
	}

	uint32_t fwd;
	uint32_t rev;
};

kmer_t rev_comp(kmer_t in, uint32_t size);
kmer_t canonicalize(kmer_t in, int size);
kmer_t canonicalize(kmer_t in, int size, bool& flipped);

inline kmer_t make_kmer(const dna_const_iterator& it, int size)
{
	kmer_t out = 0;
	for (int i = 0; i < size; i++) {
		out <<= 2;
		out |= (int) *(it + i);
	}
	return out;			
}

inline kmer_t left(const kmer_t& x, int size, int left_size)
{ 
	return x >> ((size - left_size) * 2); 
}

inline kmer_t right(const kmer_t& x, int right_size)
{ 
	return x & ((1ULL << (2 * right_size)) - 1ULL); 
}

inline kmer_t append(const kmer_t& k1, const kmer_t& k2, int size)
{ 
	return (k1 << (2 * size)) | k2; 
}

inline void rotate_left(kmer_t& x, int size, int& base)
{ 
	int base_out = left(x, size, 1); 
	x = append(right(x, size - 1), base, 1); 
	base = base_out; 
}

inline void rotate_right(kmer_t& x, int size, int& base)
{ 
	int base_out = right(x, 1); 
	x = append(base, left(x, size, size - 1), size - 1); 
	base = base_out;
}

inline
unsigned kmer_bit_value(char ascii) {
	switch (ascii) {
		case 'a':
		case 'A':
			return 0;
		case 'c':
		case 'C':
			return 1;
		case 'g':
		case 'G':
			return 2;
		case 't':
		case 'T':
			return 3;
		default:
			return 0;
	}
}

class kmer_str_iterator : public std::iterator<std::forward_iterator_tag, kmer_t>
{
public:
	explicit kmer_str_iterator(const std::string& seq, size_t kmer_size, size_t offset, kmer_t kmer)
		: m_seq(seq)
		, m_offset(offset)
		, m_kmer_size(kmer_size)
		, m_kmer(kmer)
	{}

	kmer_t operator *() { return m_kmer; }

	kmer_str_iterator& operator ++()
	{
		m_kmer <<= 2;
		m_kmer &= (1ULL << (2ULL * m_kmer_size)) - 1ULL;
		m_kmer |= kmer_bit_value(m_seq[m_offset]);
		m_offset++;
		return *this;
	}

	bool operator ==(const kmer_str_iterator& rhs) const { return m_offset == rhs.m_offset; }
	bool operator !=(const kmer_str_iterator& rhs) const { return m_offset != rhs.m_offset; }

private:
	const std::string& m_seq;
	size_t m_offset;
	const size_t m_kmer_size;
	kmer_t m_kmer;
};

// kmer_str_view presents a view of all of the kmers contained within
// a DNA seuence represented by a string.  It will have seq.size() + 1
// - kmer_size elements.
class kmer_str_view
{
public:
	kmer_str_view(const std::string& seq, size_t kmer_size)
		: m_seq(seq) 
		, m_kmer_size(kmer_size)
	{}

	kmer_str_iterator begin() const
	{
		kmer_t kmer = 0;
		for (size_t i = 0; i < m_kmer_size; i++) {
			kmer <<= 2;
			kmer |= kmer_bit_value(m_seq[i]);
		}
		return kmer_str_iterator(m_seq, m_kmer_size, m_kmer_size, kmer);
	}

	kmer_str_iterator end() const
	{
		return kmer_str_iterator("", 0, m_seq.size() + 1, 0);
	}

private:
	const std::string& m_seq;
	const size_t m_kmer_size;
};

inline kmer_t operator"" _kmer(const char* str, size_t sz) {
	return dna_sequence(std::string(str, sz)).as_kmer();
}

inline std::string kmer_str(kmer_t k, size_t s) {
	return dna_sequence(k, s).as_string();
}

class kmer_view {
 public:
  class iterator {
   public:
    iterator(dna_const_iterator cur_pos, kmer_t cur_kmer, size_t kmer_size)
        : m_cur_pos(cur_pos), m_cur_kmer(cur_kmer), m_kmer_size(kmer_size) {}
    kmer_t operator*() const {
      return right((m_cur_kmer << 2) | int(*m_cur_pos), m_kmer_size);
    }
    iterator& operator++() {
      m_cur_kmer <<= 2;
      m_cur_kmer |= int(*m_cur_pos);
      ++m_cur_pos;
      return *this;
    }

    bool operator==(const iterator& rhs) const {
      return m_cur_pos == rhs.m_cur_pos;
    }
    bool operator!=(const iterator& rhs) const {
      return m_cur_pos != rhs.m_cur_pos;
    }
    bool operator-(const iterator& rhs) const {
      return m_cur_pos - rhs.m_cur_pos;
    }

   private:
    dna_const_iterator m_cur_pos;
    kmer_t m_cur_kmer = 0;
    unsigned m_kmer_size = 0;
  };

  kmer_view(dna_slice seq, size_t kmer_size)
      : m_seq(seq), m_kmer_size(kmer_size) {
    CHECK_GT(kmer_size, 0);
  }
  kmer_view() = delete;

  iterator begin() const {
    if (m_seq.size() < m_kmer_size) {
      return end();
    }

    iterator first(m_seq.begin(), int(*m_seq.begin()), m_kmer_size);
    for (size_t i = 0; i < (m_kmer_size - 1); ++i) {
      ++first;
    }
    return first;
  }

  iterator end() const { return iterator(m_seq.end(), 0, m_kmer_size); }

 private:
  dna_slice m_seq;
  size_t m_kmer_size = 0;
};

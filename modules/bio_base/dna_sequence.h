#pragma once

#include "base/base.h"
#include "modules/bio_base/dna_base.h"
#include "modules/io/transfer_object.h"
#include "modules/io/log.h"

#include <vector>
#include <string>
#include <iostream>
#include <boost/iterator/iterator_facade.hpp>

typedef uint64_t kmer_t;

template<bool is_const>
struct type_decode;
				
// A proxy object to a binary encoded dna_base
class dna_base_proxy 
	: public boost::totally_ordered<dna_base> // Impements >, <=, !=, etc from < and ==
{
	friend struct type_decode<false>;
	friend class dna_slice;
	dna_base_proxy();  // No default constructor
	dna_base_proxy(unsigned char* data, int shift, bool rev_comp)
		: m_data(data), m_shift(shift), m_rev_comp(rev_comp) {}
public:
	operator dna_base() const 
	{
		int base = m_rev_comp ? (3-get()) : get();
		return dna_base(base);
	}
	explicit operator int() const
	{
		return m_rev_comp ? (3-get()) : get();
	}
	explicit operator char() const
	{
		return "ACGT"[m_rev_comp ? (3-get()) : get()];
	}
	void operator=(const dna_base& x) const 
	{
		int base = m_rev_comp ? (3 - int(x)) : int(x);
		put(base);	
	}
	void operator=(const dna_base_proxy& rhs) const
	{
		put(rhs.get());
	}

	dna_base complement() const { return dna_base(*this).complement(); }

private:
	unsigned char* m_data;
	int m_shift;
	bool m_rev_comp;
	int get() const { return ((*m_data) >> m_shift) & 3; }
	void put(int base) const { *m_data &= ~(3 << m_shift); *m_data |= base << m_shift; }
};

inline bool operator==(const dna_base_proxy& a, const dna_base_proxy& b) { return int(a) == int(b); }
inline bool operator<(const dna_base_proxy& a, const dna_base_proxy& b) { return int(a) < int(b); }

// The result of a comparison between two DNA sequences.
enum class dna_compare_result {
  FIRST_IS_LESS,
  FIRST_IS_PREFIX,

  EQUAL,

  SECOND_IS_PREFIX,
  SECOND_IS_LESS
};

std::ostream& operator<<(std::ostream& os, dna_compare_result compare_result);

template<bool is_const>
struct type_decode
{
	typedef unsigned char* data_type;
	typedef dna_base value_type;
	typedef dna_base_proxy ref_type;
	static ref_type do_deref(data_type data, ptrdiff_t offset, bool complement = false)
	{
		return dna_base_proxy(data + offset/4, 6-(offset & 3)*2, complement); 
	}
};

template<>
struct type_decode<true>
{
	typedef const unsigned char* data_type;
	typedef const dna_base value_type;
	typedef dna_base ref_type;
	static ref_type do_deref(data_type data, ptrdiff_t offset, bool complement = false)
	{
		int base = (data[offset/4] >> (6-(offset&3)*2)) & 3;
		return dna_base(complement ? 3 - base : base);
	}

	static dna_compare_result do_compare(data_type data1, ptrdiff_t offset1, bool rev_comp1,
                                      data_type data2, ptrdiff_t offset2, bool rev_comp2,
                                      size_t len1, size_t len2);

};

template<bool is_const>
class dna_iterator_temp : 
	public boost::iterator_facade<
		dna_iterator_temp<is_const>,
		typename type_decode<is_const>::value_type,
		std::random_access_iterator_tag,
		typename type_decode<is_const>::ref_type>
{
protected:
	typedef typename type_decode<is_const>::data_type data_type;
	typedef typename type_decode<is_const>::value_type val_type;
	typedef typename type_decode<is_const>::ref_type ref_type;

	friend class boost::iterator_core_access;	
public:
	dna_iterator_temp() 
		: m_data(NULL), m_offset(0), m_rev_comp(false) {}
	dna_iterator_temp(const dna_iterator_temp<false>& rhs) 
		: m_data(data_type(rhs.get_data())), m_offset(rhs.get_offset()), m_rev_comp(rhs.is_rev_comp()) {}
	dna_iterator_temp(data_type data, ptrdiff_t offset, bool rev_comp)
		: m_data(data), m_offset(offset), m_rev_comp(rev_comp) {}
	dna_iterator_temp rev_comp() const { return dna_iterator_temp(m_data, m_offset, !m_rev_comp); }
	data_type get_original_data() const { return m_data; }
	data_type get_data() const { return m_data; }
	ptrdiff_t get_offset() const { return m_offset; }
	ptrdiff_t get_original_offset() const { return m_offset; }
	ptrdiff_t get_offset_in_byte() const { return m_offset & 3; }
	data_type get_data_byte() const { return m_data + m_offset/4; }

	bool is_rev_comp() const { return m_rev_comp; }
	
private:
	void increment() { if (m_rev_comp) m_offset--; else m_offset++; }
	void decrement() { if (m_rev_comp) m_offset++; else m_offset--; }
	void advance(ptrdiff_t n) { if (m_rev_comp) m_offset -= n; else m_offset += n; }
	ptrdiff_t distance_to(const dna_iterator_temp& rhs) const 
	{
		if (m_rev_comp) {
			CHECK(rhs.m_rev_comp);
			return m_offset - rhs.m_offset;
		} else {
			CHECK(!rhs.m_rev_comp);
			return rhs.m_offset - m_offset; 
		}
	}
	bool equal(const dna_iterator_temp& rhs) const { return m_data == rhs.m_data && m_offset == rhs.m_offset; }	
	ref_type dereference() const 
		{ return type_decode<is_const>::do_deref(m_data, m_offset, m_rev_comp); }
	data_type m_data;
	ptrdiff_t m_offset;
	bool m_rev_comp;	
};

typedef dna_iterator_temp<false> dna_iterator;
typedef dna_iterator_temp<true> dna_const_iterator;

// No bound checks are done. At least len bases must be defined in the
// original dna_sequence instances.
inline dna_compare_result subseq_compare(dna_const_iterator start1,
                                         dna_const_iterator start2, size_t len1,
                                         size_t len2) {
  if (start1 == start2 && start1.is_rev_comp() == start2.is_rev_comp()) {
    // If we're operating on large repositories, it can be slow to dereference
    // start1 and start2.  Avoid it if we can...
    if (len1 < len2) {
      return dna_compare_result::FIRST_IS_PREFIX;
    } else if (len1 > len2) {
      return dna_compare_result::SECOND_IS_PREFIX;
    } else {
      return dna_compare_result::EQUAL;
    }
  }
  return type_decode<true>::do_compare(
      start1.get_data_byte(), start1.get_offset_in_byte(), start1.is_rev_comp(),
      start2.get_data_byte(), start2.get_offset_in_byte(), start2.is_rev_comp(),
      len1, len2);
}

inline bool subseq_equal(dna_const_iterator start1, dna_const_iterator start2,
                         size_t len) {
  switch (subseq_compare(start1, start2, len, len)) {
    case dna_compare_result::EQUAL:
      return true;
    default:
      return false;
  }
}

inline bool subseq_lessthan(dna_const_iterator start1,
                            dna_const_iterator start2, size_t len1,
                            size_t len2) {
  switch( subseq_compare(start1, start2, len1, len2)) {
    case dna_compare_result::FIRST_IS_LESS:
    case dna_compare_result::FIRST_IS_PREFIX:
      return true;
    default:
      return false;
  }
}

class dna_slice;
class dna_sequence
	: boost::totally_ordered<dna_sequence>
	, boost::addable<dna_sequence>
{
public:
	typedef dna_const_iterator const_iterator;
	typedef dna_iterator iterator;

	dna_sequence() = default;  // Construct empty sequence
	dna_sequence(const dna_sequence& rhs);  // Copy constructor
	dna_sequence& operator=(const dna_sequence& rhs);  // Assignment operator
	dna_sequence& operator=(const dna_slice& rhs);  // Assignment operator
	dna_sequence(dna_sequence&& rhs);  // Move constructor
	dna_sequence& operator=(dna_sequence&& rhs);  // Move assignment operator
	~dna_sequence(); // Destructor

	// Conversion + convience constructors
	explicit dna_sequence(const const_iterator& start, const const_iterator& end); 
	explicit dna_sequence(size_t size); // Construct a sequence of fixed size 
	explicit dna_sequence(const dna_base& b, size_t count = 1); // Construct from a base
	explicit dna_sequence(kmer_t kmer, int size);  // Construct from kmer_t
	explicit dna_sequence(const dna_slice& rhs);  // Convert from slice

	void clear();
	void push_back(dna_base);

	using size_type = size_t;
	size_t size() const { return m_size; }
	size_t isize() const { return m_size/4 + 1; }
	bool empty() const { return size() == 0; }

	// Deserialize version
	dna_sequence(const std::string& seq, bool packed = false);

	const dna_base operator[](uint64_t i) const
	{ return type_decode<true>::do_deref(m_data, i + 1); }
	dna_base_proxy operator[](uint64_t i) 
	{ return type_decode<false>::do_deref(m_data, i + 1); }
	
	iterator begin() { return iterator(m_data, 1, false); }
	const_iterator begin() const { return const_iterator(m_data, 1, false); }
	iterator end() { return iterator(m_data, 1 + m_size, false); }
	const_iterator end() const { return const_iterator(m_data, 1 + m_size, false); }

	iterator rcbegin() { return iterator(m_data, m_size, true); }
	const_iterator rcbegin() const { return const_iterator(m_data, m_size, true); }
	iterator rcend() { return iterator(m_data, 0, true); }
	const_iterator rcend() const { return const_iterator(m_data, 0, true); }
	
	void reserve(uint64_t new_capacity);

	std::string as_string() const;
	std::string as_packed() const;
	kmer_t as_kmer() const;

	bool operator<(const dna_sequence& rhs) const;
	bool operator==(const dna_sequence& rhs) const;
	bool operator==(const dna_slice& rhs) const;
	
	dna_sequence subseq(size_t offset, size_t len) const;
	dna_sequence reverse() const; // Just reverse, no complement
	dna_sequence rev_comp() const; // Reverse complement
	dna_sequence canonicalize() const { bool ignore; return canonicalize(ignore); }
	dna_sequence canonicalize(bool& flipped) const;  // RComplement or not...

	dna_sequence& operator+=(const dna_sequence& rhs);
	dna_sequence& operator+=(const dna_base& rhs);
	dna_sequence& operator+=(const dna_slice& rhs);

	dna_compare_result compare_to(const dna_slice& rhs) const;
	unsigned shared_prefix_length(const dna_slice& rhs) const;

	static iterator copy_bases(const dna_slice& slice, iterator dest);

private:
	void resize(size_t new_size);
	unsigned char* m_data = nullptr; // Bit encoded data
	size_t m_capacity = 0; // Capacity in bases
	size_t m_size = 0;   // Size in bases
};

class dna_slice
	: boost::totally_ordered<dna_slice>
{
public:
	typedef dna_const_iterator const_iterator;
	dna_slice() : m_size(0) {} // Make empty slice
	dna_slice(const dna_sequence& rhs) 
		: m_begin(rhs.begin()), m_size(rhs.size()) {}  // Seq->Slice
	dna_slice(const const_iterator& begin, size_t size)
		: m_begin(begin), m_size(size) {}  // Generate from iterator + size
	dna_slice(const const_iterator& begin, const const_iterator& end)
		: m_begin(begin), m_size(end - begin) {}  // Generate from iterators
	// Copy constructor and assignment are default

	size_t size() const { return m_size; }
	bool empty() const { return size() == 0; }

	dna_base operator[](uint64_t i) const { return *(m_begin + i); }
	dna_base_proxy operator[](uint64_t i) { 
		dna_const_iterator it = m_begin + i;
		unsigned char* data = const_cast<unsigned char*>(it.get_data());
		size_t offset = it.get_offset();
		return dna_base_proxy(data + offset/4, 6-(offset & 3)*2, it.is_rev_comp());
	}
	
	const_iterator begin() const { return m_begin; }
	const_iterator end() const { return m_begin + m_size; }

	const_iterator rcbegin() const { return (m_begin + (m_size - 1)).rev_comp(); }
	const_iterator rcend() const { return (m_begin - 1).rev_comp(); }

	std::string as_string() const;
	kmer_t as_kmer() const;

	bool operator<(const dna_slice& rhs) const {
		return subseq_lessthan(m_begin, rhs.m_begin, m_size, rhs.m_size); 
	}
	bool operator==(const dna_slice& rhs) const
	{
		if (m_size != rhs.m_size) return false;
		return subseq_equal(m_begin, rhs.m_begin, m_size); 
	}

	dna_compare_result compare_to(const dna_slice& rhs) const {
      return subseq_compare(m_begin, rhs.m_begin, m_size, rhs.m_size);
    }
	unsigned shared_prefix_length(const dna_slice& rhs) const;
	
	dna_slice subseq(size_t offset, size_t len) const {
		CHECK_LE(offset + len, m_size) << "Offset: " << offset << " Len: " << len;
		return dna_slice(m_begin + offset, len);
	}
	dna_slice rev_comp() const { return dna_slice(rcbegin(), rcend()); }
	dna_slice canonicalize() const { bool ignore; return canonicalize(ignore); }
	dna_slice canonicalize(bool& flipped) const {
		if (*this < this->rev_comp()) { flipped = false; return *this; }
		else { flipped = true; return this->rev_comp(); }
	} 

private:
	const_iterator m_begin;
	size_t m_size;
};

inline dna_compare_result dna_sequence::compare_to(const dna_slice& rhs) const {
  return dna_slice(*this).compare_to(rhs);
}

inline unsigned dna_sequence::shared_prefix_length(const dna_slice& rhs) const {
  return dna_slice(*this).shared_prefix_length(rhs);
}

// Make dna sequences easily printable in error messages, etc.
std::ostream& operator<<(std::ostream&, const dna_sequence&);
std::ostream& operator<<(std::ostream&, const dna_slice&);

// Function used to print dna sequences.  Tests can modify this to get
// additional debugging information outputted.
namespace dna_testutil {
std::ostream& default_dna_printer(std::ostream&,
                                  const dna_slice&);
extern std::function<std::ostream&(std::ostream&, const dna_slice&)>
    g_dna_printer;
}  // namespace dna_testutil

template<>
struct transfer_info<dna_sequence>
{
	typedef std::string type;
	static std::string get(const dna_sequence& value)
	{
		return value.as_packed(); 
	}
	static void put(dna_sequence& value, const std::string& str)
	{
		value = dna_sequence(str, true);
	}
};

uint8_t byte_rev_comp_bases(uint8_t a_byte); // complement an entire byte.
uint64_t long_rev_comp_bases(uint64_t a_byte); // complement a 64 bit long.

// Helper class for DNA ordering; supplies < and == operators based on
// a "compare_to" call that returns a dna_compare_result.
//
// TODO(nils): Can we have dna_sequence and dna_slice use this?
template <typename this_type, typename target_type>
class dna_sequence_ordered {
 public:
  friend bool operator==(const this_type& lhs, const target_type& rhs) {
    dna_compare_result cmp = lhs.compare_to(rhs);
    return cmp == dna_compare_result::EQUAL;
  }
  friend bool operator<(const this_type& lhs, const target_type& rhs) {
    dna_compare_result cmp = lhs.compare_to(rhs);
    return cmp == dna_compare_result::FIRST_IS_LESS ||
           cmp == dna_compare_result::FIRST_IS_PREFIX;
  }
  friend bool operator>(const this_type& lhs, const target_type& rhs) {
    dna_compare_result cmp = lhs.compare_to(rhs);
    return cmp == dna_compare_result::SECOND_IS_LESS ||
           cmp == dna_compare_result::SECOND_IS_PREFIX;
  }
  friend bool operator<=(const this_type& lhs, const target_type& rhs) {
    return !(lhs > rhs);
  }
  friend bool operator>=(const this_type& lhs, const target_type& rhs) {
    return !(lhs < rhs);
  }
  friend bool operator!=(const this_type& lhs, const target_type& rhs) {
    return !(lhs == rhs);
  }

  bool is_prefix_or_equal(target_type& rhs) const {
    dna_compare_result cmp =
        static_cast<const this_type&>(*this).compare_to(rhs);
    return cmp == dna_compare_result::FIRST_IS_PREFIX ||
           cmp == dna_compare_result::EQUAL;
  }
};

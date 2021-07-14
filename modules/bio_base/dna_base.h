
#ifndef __dna_base_h__
#define __dna_base_h__

#include <boost/operators.hpp>
#include <ostream>

#include "modules/io/io.h"
#include "modules/io/utils.h"

namespace internal {
class dna_all_bases_t;
class dna_all_bases_iterator;
}  // namespace internal

// Make a dna_base a first class type
class dna_base 
	: public boost::totally_ordered<dna_base> // Impements >, <=, !=, etc from < and ==
{
public:
	dna_base() : m_base(0) {}
	explicit dna_base(char b);
	explicit dna_base(int b);

	explicit operator char() const { return "ACGT"[m_base]; }
	explicit operator int() const { return m_base; }

	dna_base complement() const { return dna_base(3-m_base); } // A<->T, C<->G

	static constexpr int k_num_bases = 4;

private:
	constexpr dna_base(int b, bool) : m_base(b) {}
	friend class internal::dna_all_bases_iterator;
	unsigned char m_base;
};

inline dna_base::dna_base(char c)
{
	switch(c)
	{
	case 'a':
	case 'A':
		m_base = 0; break;
	case 'c':
	case 'C':
		m_base = 1; break;
	case 'g':
	case 'G':
		m_base = 2; break;
	case 't':
	case 'T':
		m_base = 3; break;
	default:
		throw io_exception(printstring("Failed conversion of dna_base, c = '%c'", c));
	}
}

namespace internal {

class dna_all_bases_iterator {
 public:
  dna_all_bases_iterator operator++() {
    m_b++;
    return *this;
  }
  bool operator!=(const dna_all_bases_iterator &rhs) const {
    return m_b != rhs.m_b;
  }
  dna_base operator*() const { return dna_base(m_b, true); }

 private:
  friend class dna_all_bases_t;
  constexpr dna_all_bases_iterator(int b) : m_b(b) {}
  int m_b;
};

// Implement the bare minimum required so we can use this in a range
// for loop.
class dna_all_bases_t {
 public:
  constexpr dna_all_bases_iterator begin() const {
    return dna_all_bases_iterator(0);
  }
  constexpr dna_all_bases_iterator end() const  {
    return dna_all_bases_iterator(4);
  }
};

}  // namespace internal

// Allow easy iteration through all bases by:
// for (dna_base b : dna_bases()) {
//   ...
// }
inline constexpr internal::dna_all_bases_t dna_bases() {
  return internal::dna_all_bases_t();
}

// This is a fixed size array that uses a dna_base as the array index.
template <typename T>
class dna_base_array : public std::array<T, dna_base::k_num_bases> {
  using base_array = std::array<T, dna_base::k_num_bases>;
 public:
  T& operator[](dna_base b) { return (*(base_array*)this)[int(b)]; }
  const T& operator[](dna_base b) const {
    return (*(const base_array*)this)[int(b)];
  }
};


inline bool operator==(const dna_base& a, const dna_base& b) { return int(a) == int(b); }
inline bool operator<(const dna_base& a, const dna_base& b) { return int(a) < int(b); }

// Make a dna_base a first class type
class dna_del_base 
	: public boost::totally_ordered<dna_del_base> // Impements >, <=, !=, etc from < and ==
{
public:
	dna_del_base() : m_base(0) {}
	explicit dna_del_base(char b);
	explicit dna_del_base(int b);

	explicit operator char() const { return "ACGT."[m_base]; }
	explicit operator int() const { return m_base; }

	bool operator==(const dna_del_base& rhs) const { return m_base == rhs.m_base; }
	bool operator<(const dna_del_base& rhs) const { return m_base < rhs.m_base; }

	dna_del_base complement() const { 
		if (m_base == 4) return *this;
		return dna_del_base(3-m_base); 
	} // A<->T, C<->G

private:
	unsigned char m_base;
};

inline dna_base::dna_base(int b)
{
	if (b < 0 || b > 3)
		throw io_exception("Conversion from int to dna_base failed");
	m_base = b;
}

inline std::ostream& operator<<(std::ostream& os, dna_base b) {
  return os << char(b);
}

#endif


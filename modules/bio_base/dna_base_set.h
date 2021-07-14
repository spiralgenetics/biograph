
#ifndef __dna_base_set__
#define __dna_base_set__

#include <boost/operators.hpp>
#include "modules/bio_base/dna_base.h"

class dna_base_set
	: public boost::andable<dna_base_set>  // provides a & b from a &= b
	, public boost::orable<dna_base_set>
	, public boost::subtractable<dna_base_set>
	, public boost::totally_ordered<dna_base_set>
{
public:
	dna_base_set() : m_set(0) {}  // Empty set (ie '-')
	dna_base_set(char c);  // User IUPAC Nucleic Acid Code
	dna_base_set(const dna_base& base) : m_set(1 << ((int) base)) {}
	dna_base_set& operator|=(const dna_base_set& rhs) 
		{ m_set |= rhs.m_set; return *this; } // Union
	dna_base_set& operator&=(const dna_base_set& rhs) 
		{ m_set &= rhs.m_set; return *this; } // Intersection
	dna_base_set& operator-=(const dna_base_set& rhs) 
		{ m_set &= ~(rhs.m_set); return *this; } // Set difference
	bool operator<(const dna_base_set& rhs) const { return m_set < rhs.m_set; }
	bool operator==(const dna_base_set& rhs) const { return m_set == rhs.m_set; }
	operator bool() { return m_set != 0; } // True if non-empty

// reduce fn over the set s, with 'result' initially set to the first, potentially empty value.
// Implies that fn is a binary function or functor with the following signature: void fn( const dna_base&, T& )
	template <typename T, typename Functor> static void reduce(const Functor& fn, const dna_base_set& s, T& result)
		{
			if (s & dna_base_set('A')) fn(dna_base('A'), result);
			if (s & dna_base_set('C')) fn(dna_base('C'), result);
			if (s & dna_base_set('G')) fn(dna_base('G'), result);
			if (s & dna_base_set('T')) fn(dna_base('T'), result);
		}
	std::string as_list(char sep = ',');  // sep = 0 means no seperator
	char as_code();  // Uses IUPAC Nucleic Acid Code
private:
	unsigned int m_set;  // Bit encoded set
};

// Reverses the passed string in place.  Copy the input if you need to keep it.
void reverse_complement_iupac_string(std::string& iupac_string);

#endif

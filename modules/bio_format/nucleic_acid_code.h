#ifndef __nucleic_acid_code_H__
#define __nucleic_acid_code_H__

#include <map>
#include <algorithm>
#include <string.h>
// IUPAC Nucleic Acid Code
// as described here: http://www.dna.affrc.go.jp/misc/MPsrch/InfoIUPAC.html
// Usage: NAC['B'] -> "GTCU"

class IUPAC_NAC : public std::map<char,const char*>
{
public:
	static IUPAC_NAC& instance() { static IUPAC_NAC _nac; return _nac; }
private:
	IUPAC_NAC()
		{
			(*this)['A'] = "A";
			(*this)['T'] = "T";
			(*this)['G'] = "G";
			(*this)['C'] = "C";
			(*this)['U'] = "U";
			(*this)['R'] = "GA";
			(*this)['Y'] = "TCU";
			(*this)['K'] = "GTU";
			(*this)['M'] = "AC";
			(*this)['S'] = "GC";
			(*this)['W'] = "ATU";
			(*this)['B'] = "GTCU";
			(*this)['D'] = "GATU";
			(*this)['H'] = "ATCU";
			(*this)['V'] = "GCA";
			(*this)['N'] = "ATGCU";
			(*this)['-'] = "";
		};
	IUPAC_NAC( const IUPAC_NAC& );
	IUPAC_NAC& operator=( const IUPAC_NAC& );
};

#define NAC IUPAC_NAC::instance()

#endif //__nucleic_acid_code_h__


#include "modules/bio_base/dna_base.h"
#include "modules/io/utils.h"

dna_del_base::dna_del_base(char c)
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
	case '.':
		m_base = 4; break;
	default:
		throw io_exception(printstring("Failed conversion of dna_base, c = '%c'", c));
	}
}

dna_del_base::dna_del_base(int b)
{
	if (b < 0 || b > 4)
		throw io_exception("Conversion from int to dna_base failed");
	m_base = b;
}



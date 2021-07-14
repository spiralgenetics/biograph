
#ifndef __dna_multiseq_h__
#define __dna_multiseq_h__

#include "modules/bio_base/dna_sequence.h"

class dna_multiseq
{
public:
	dna_multiseq(const dna_sequence& s1, const dna_sequence& s2);
	std::string get_string(size_t which);
	
	//std::string as_string();
	//std::string as_packed();
private:
	typedef std::vector<dna_del_base> dna_del_seq;
	typedef std::vector<dna_del_seq> del_seqs_t;
	del_seqs_t m_seqs;
	
};

#endif

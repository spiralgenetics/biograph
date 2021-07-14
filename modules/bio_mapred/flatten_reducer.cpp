
#include "modules/bio_mapred/flatten_reducer.h"

REGISTER_1(reducer, flatten, const std::string&);

void flatten_reducer::typed_start(const dna_sequence& key)
{
	if (m_cur.size() == 0) return;
	int shared;
	flatten_value fv;
	make_output(shared, fv, m_last.size() == 0);
	output(shared, fv);
}

void flatten_reducer::typed_add_value(const dna_sequence& key, int value) 
{ 
	if (value == -1) {
		m_cur = key;
	}
	else {
		m_bits |= (1 << value);
	} 
}

void flatten_reducer::typed_end() 
{
}

void flatten_reducer::finalize(kv_sink& context)
{
	int shared;
	flatten_value fv;
	make_output(shared, fv, m_last.size() == 0);
	context.write_msgpack(shared, fv);
}

void flatten_reducer::make_output(int& shared, flatten_value& fv, bool include_seq)
{
	if (m_last.size() == 0) {
		shared = -1;
	} else {
		shared = 0;
		for(size_t i = 0; i < std::min(m_last.size(), m_cur.size()); i++)
		{
			if (m_last[i] != m_cur[i]) break;
			shared++;
		}
	}
	fv.context = m_cur.size();
	fv.bits = m_bits;
	if (include_seq)
		fv.seq = m_cur;

	/*
	static int ii = 0;
	std::string line = printstring("%d: %s: ", ii++, m_cur.as_string().c_str());
	for(int i = 0; i < 4; i++)
		if (m_bits & (1<<i)) line += printstring("%c ", (char) dna_base(i));
	SPLOG("%s", line.c_str());
	*/

	m_last = m_cur;
	m_bits = 0;
}

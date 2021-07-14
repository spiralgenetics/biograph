
#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/mapred/reducer.h"

// key is: int shared;  // How much do I share with next seq, -1 for first sequence in file

struct flatten_value
{
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(seq);
		FIELD(context);
		FIELD(bits);
	}
	dna_sequence seq;   // Only set of first and last record of each manifest
	uint8_t context;  // How long is the entry
	uint8_t bits; // Bit field 1,2,4,8 = A,C,G,T
};

class flatten_reducer : public typed_reducer<flatten_reducer, dna_sequence, int, int, flatten_value>
{
public:
	flatten_reducer(const std::string& params) : m_bits(0) {}

	void typed_start(const dna_sequence& key);
	void typed_add_value(const dna_sequence& key, int value);
	void typed_end();
	void finalize(kv_sink& context) override;

private:
	void make_output(int& shared, flatten_value& fv, bool include_seq);
	dna_sequence m_last;
	dna_sequence m_cur;
	uint8_t m_bits;
};


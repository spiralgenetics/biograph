#pragma once

#include <bitset>
#include <deque>
#include "modules/io/prefix_sum.h"
#include "modules/io/range_coder.h"

class dict_codec
{
public:
	// Number of bits in dictionary, must be  8 < x < 16
	dict_codec(size_t dict_bits);
	// Zero out dictionary
	void reset();

	// Encode + decode.  If you need to write with record boundries,
	// you must call eor() during encode and at the same times on decode
	// eor() makes sure that all data needed to decode a record is available
	// and not buffered in the dict_codec
	void encode(range_encoder& r, uint8_t byte);
	void enc_eor(range_encoder& r);
	void enc_eof(range_encoder& r);

	bool decode(range_decoder& r, uint8_t& byte);
	void dec_eor(range_decoder& r);
private:
	static const uint16_t k_nil_entry = 0;
	static const uint32_t k_dist_mult = 50;

	size_t m_dict_bits;
	size_t m_dict_size;

	// Current dictionary size
	size_t m_cur_size;

	// Current entry
	uint16_t m_cur_entry;

	// Delay follows info
	uint16_t m_old_entry;
	uint8_t m_old_byte;

	// Actual dictionary entries
	uint16_t get_parent(uint16_t entry);
	uint8_t get_byte(uint8_t entry);
	std::vector<uint16_t> m_parent;
	std::vector<uint8_t> m_byte;

	// Hash based lookup table for finding entries
	uint16_t hash_lookup(uint16_t parent, uint8_t byte);
	// Add/lookup an entry, returns true if added, false is already an entry
	bool hash_put(uint16_t parent, uint8_t byte, uint16_t entry);
	std::vector<uint32_t> m_hash_key;
	std::vector<uint16_t> m_hash_value;

	// Tracking for new potential entries
	void add_follows(uint16_t entry, uint8_t byte);
	typedef std::bitset<256> bits_t;
	std::vector<bits_t> m_follows;

	// Method to compute symbols probabilities
	prefix_sum_dist m_dist;

	// Used for decoding state
	std::deque<uint8_t> m_decode_buf;
};

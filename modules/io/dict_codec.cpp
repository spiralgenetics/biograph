#include "modules/io/dict_codec.h"

// We preallocate everything since we know max sizes in advance, and they are small
dict_codec::dict_codec(size_t dict_bits)
	: m_dict_bits(dict_bits)
	, m_dict_size(1 << dict_bits)
	, m_cur_size(257)
	, m_cur_entry(0)
	, m_old_entry(0)
	, m_old_byte(0)
	, m_parent(m_dict_size)
	, m_byte(m_dict_size)
	, m_hash_key(4*m_dict_size)
	, m_hash_value(4*m_dict_size)
	, m_follows(m_dict_size)
	, m_dist(m_dict_size)
{
	// Zap the hashes to invalid keys
	for(size_t i = 0; i < 4*m_dict_size; i++) {
		m_hash_key[i] = 0xffffffff;
	}
	// Fill in an initial dictionary of all values
	for(size_t i = 0; i < 256; i++) {
		m_byte[i+1] = uint8_t(i);
		hash_put(0, uint8_t(i), i+1);
		m_dist.inner().add(i+1, 1);
	}
	m_dist.inner().add(0, 1);
}

void dict_codec::reset()
{
	// SPLOG("Doing a reset");
	m_cur_size = 257;
	m_cur_entry = 0;
	m_old_entry = 0;
	m_old_byte = 0;
	m_dist.inner().reset();
	for(size_t i = 0; i < m_dict_size; i++) {
		m_follows[i].reset();
	}
	for(size_t i = 0; i < 4*m_dict_size; i++) {
		m_hash_key[i] = 0xffffffff;
		m_hash_value[i] = 0;
	}
	for(size_t i = 0; i < 256; i++) {
		hash_put(0, uint8_t(i), i+1);
		m_dist.inner().add(i+1, 1);
	}
	m_dist.inner().add(0, 1);
}

void dict_codec::encode(range_encoder& r, uint8_t byte)
{
	// SPLOG("Dict = %d, byte = %d", m_cur_entry, byte);
	// Check if we can 'continue' the current entry
	uint16_t next = hash_lookup(m_cur_entry, byte);
	if (next == k_nil_entry) {
		// Nope, Write current entry out
		// SPLOG("Transfer %d", m_cur_entry);
		r.encode(m_dist, m_cur_entry);
		// Make it more likely next time
		// SPLOG("Adding %d", m_cur_entry);
		m_dist.inner().add(m_cur_entry, k_dist_mult);
		// Maybe build a new entry
		if (m_old_entry) {
			add_follows(m_old_entry, m_old_byte);
		}
		// Setup delay
		m_old_entry = m_cur_entry;
		m_old_byte = byte;
		// Add next entry
		m_cur_entry = hash_lookup(k_nil_entry, byte);
	} else {
		// Add to existing
		m_cur_entry = next;
	}
}

void dict_codec::enc_eor(range_encoder& r)
{
	if (m_cur_entry == k_nil_entry) return;
	// SPLOG("Transfer %d", m_cur_entry);
	r.encode(m_dist, m_cur_entry);
	// Make it more likely next time
	// SPLOG("Adding %d", m_cur_entry);
	m_dist.inner().add(m_cur_entry, k_dist_mult);
	// SPLOG("End of record");
	// Zero our state
	m_old_entry = k_nil_entry;
	m_cur_entry = k_nil_entry;
	if (m_dist.inner().total(m_dict_size) > k_min_range/2) {
		reset();
	}
}

void dict_codec::enc_eof(range_encoder& r)
{
	// flush it
	enc_eor(r);

	r.encode(m_dist, 0);
}

bool dict_codec::decode(range_decoder& r, uint8_t& byte)
{
	// If we don't have any buffered data to dump
	if (m_decode_buf.size() == 0) {
		if (m_old_entry != 0) {
			add_follows(m_old_entry, m_old_byte);
		}
		// Pull a new value
		m_old_entry = m_cur_entry;
		m_cur_entry = r.decode(m_dist);
		// SPLOG("Transfer %d", m_cur_entry);
		if(m_cur_entry == 0) {
			return false;
		}
		// Update count
		// SPLOG("Adding %d", m_cur_entry);
		m_dist.inner().add(m_cur_entry, k_dist_mult);
		// Decode into output buffer
		uint16_t it = m_cur_entry;
		while(it != k_nil_entry) {
			// SPLOG("  Pushing front: %c", m_byte[it]);
			m_decode_buf.push_front(m_byte[it]);
			it = m_parent[it];
		}
		m_old_byte = m_decode_buf.front();
	}
	byte = m_decode_buf.front();
	m_decode_buf.pop_front();
	return true;
}

void dict_codec::dec_eor(range_decoder& r)
{
	// SPLOG("End of record");
	m_cur_entry = 0;
	m_old_entry = 0;
	if (m_dist.inner().total(m_dict_size) > k_min_range/2) {
		reset();
	}
}

// Why our own hash table thing?
// 1) We never rehash, fixed maximum size
// 2) We also need 'reverse' mapping
// 2) Max performance, min memory, code transparency
// 3) Due to a bug gcc's unordered_map has perf issues in some versions
uint16_t dict_codec::hash_lookup(uint16_t parent, uint8_t byte)
{
	uint32_t key = (parent << 8) + byte;
	uint32_t hash = (key * 0x7b512cf1 + 0xbd87ad01) >> (32 - (m_dict_bits + 2));
	while(m_hash_key[hash] != key && m_hash_value[hash] != 0) {
		hash++; hash &= (4 * m_dict_size - 1);
	}
	return m_hash_value[hash];
}

bool dict_codec::hash_put(uint16_t parent, uint8_t byte, uint16_t entry)
{
	uint32_t key = (parent << 8) + byte;
	uint32_t hash = (key * 0x7b512cf1 + 0xbd87ad01) >> (32 - (m_dict_bits + 2));
	while(m_hash_key[hash] != key && m_hash_value[hash] != 0) {
		hash++; hash &= (4 * m_dict_size - 1);
	}
	if (m_hash_key[hash] == key) {
		return false;
	}
	m_hash_key[hash] = key;
	m_hash_value[hash] = entry;
	return true;
}

void dict_codec::add_follows(uint16_t entry, uint8_t byte)
{
	// If it's too big, reset
	if (m_cur_size >= m_dict_size || m_dist.inner().total(m_dict_size) > k_min_range/2) {
		reset();
		return;
	}
	if (m_follows[entry].test(byte)) {
		if (!hash_put(entry, byte, m_cur_size)) {
			// Whoops, already added
			return;
		}
		// SPLOG("New dict %d:%d -> %d", entry, byte, (int) m_cur_size);
		// Add hash + regular entry
		m_parent[m_cur_size] = entry;
		m_byte[m_cur_size] = byte;
		// Move the two prior times into new entry
		m_dist.inner().sub(entry, 2*k_dist_mult);
		m_dist.inner().add(m_cur_size, 2*k_dist_mult);
		m_cur_size++;
	} else {
		// SPLOG("Setting %d:%d", entry, byte);
		m_follows[entry].set(byte);
	}
}


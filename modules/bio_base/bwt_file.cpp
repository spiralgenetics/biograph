#include <endian.h>

#include "modules/io/msgpack_transfer.h"
#include "modules/io/make_unique.h"
#include "modules/io/version.h"
#include "modules/io/bitcount.h"
#include "modules/bio_base/bwt_file.h"

constexpr uint64_t bwt_file::k_magic;

uint64_t bwt_file::validate_offset(const char* buf)
{
	if (*reinterpret_cast<const uint64_t*>(buf) != bwt_file::k_magic) {
			throw io_exception(m_buf.path() + " does not appear to be a valid BWT.");
	}
	uint64_t offset = *reinterpret_cast<const uint64_t*>(buf + sizeof(uint64_t));
	return offset;
}

bwt_file::bwt_file(const std::string& path)
	: m_buf(path, mmap_buffer::mode::read_populate)
	, m_offset(validate_offset(m_buf.buffer()))
	, m_header(msgpack_deserialize<bwt_header>(std::string(m_buf.buffer() + m_offset, m_buf.size() - m_offset)))
	, m_base(m_buf.buffer() + (2 * sizeof(uint64_t)))
	, m_entries(m_header.ca_table[4])
	, m_bits_size(bitcount::compute_size(m_entries))
	, m_base_bits{
		{ m_base + 0 * m_bits_size, m_entries },
		{ m_base + 1 * m_bits_size, m_entries },
		{ m_base + 2 * m_bits_size, m_entries },
		{ m_base + 3 * m_bits_size, m_entries },
	}
	, m_century_bits(m_base + 4 * m_bits_size, m_entries)
	, m_century_table(reinterpret_cast<const uint32_t *>(m_base + (5 * m_bits_size)))
{
	SPLOG("bwt_file::bwt_file> BWT loaded, %lu entries", m_entries)

	// Dump all tables and bitcounts

	// printf(" # c  C A C G T\n");
	// for(size_t i = 0; i < m_entries; i++) {
	// 	printf("%2lu %d %2lu %lu %lu %lu %lu\n",
	// 		i, m_century_bits.get(i), m_century_bits.count(i),
	// 		m_base_bits[0].count(i), m_base_bits[1].count(i), m_base_bits[2].count(i), m_base_bits[3].count(i));
	// }

	// printf("\nCentury table:\n");
	// for(size_t i = 0; i < m_header.ca_table[5]; i++) {
	// 	printf("%2lu %u\n", i, m_century_table[i]);
	// }

	// printf("\nC(a) table:\n");
	// for(size_t i = 0; i < 6; i++) {
	// 	printf("%lu %lu\n", i, m_header.ca_table[i]);
	// }
}

bwt_range bwt_file::bwt() { return bwt_range(*this, 0, m_entries); }

bwt_range bwt_range::find(const std::string& query) const {
	return find(dna_sequence(query));
}

bwt_range bwt_range::find(const dna_sequence& s) const {
	size_t cur_begin = m_begin;
	size_t cur_end = m_end;
	for(int i = s.size() - 1; i >= 0; i--) {
		const dna_base& b = s[i];
		// base_start == Ca(b)
		size_t base_start = m_file.m_header.ca_table[(int) b];
		// off_begin = O(b, r)
		size_t off_begin = m_file.m_base_bits[(int) b].count(cur_begin);
		// off_end = O(b, R)
		size_t off_end = m_file.m_base_bits[(int) b].count(cur_end);
		cur_begin = base_start + off_begin;
		cur_end = base_start + off_end;
		if (cur_begin == cur_end) { break; }
	}
	return bwt_range(m_file, cur_begin, cur_end);
}

bwt_range bwt_range::push_front(const dna_base& b) const {
	// base_start == Ca(b)
	size_t base_start = m_file.m_header.ca_table[(int) b];
	// off_begin = O(b, r)
	size_t off_begin = m_file.m_base_bits[(int) b].count(m_begin);
	// off_end = O(b, R)
	size_t off_end = m_file.m_base_bits[(int) b].count(m_end);
	return bwt_range(m_file, base_start + off_begin, base_start + off_end);
}

uint32_t bwt_range::get_match(size_t which) const {
	size_t e = m_begin + which;
	// SPLOG("e == %lu", e);
	size_t dist = 0;
	while(!m_file.m_century_bits.get(e))
	{
		int b =
			(1 * m_file.m_base_bits[1].get(e)) +
			(2 * m_file.m_base_bits[2].get(e)) +
			(3 * m_file.m_base_bits[3].get(e));
		// base_start == Ca(b)
		size_t base_start = m_file.m_header.ca_table[(int) b];
		// off_begin = O(b, e)
		size_t offset = m_file.m_base_bits[(int) b].count(e);
		// compute new location
		e = base_start + offset;
		dist++;
		// SPLOG("e == %lu, dist == %lu", e, dist);
	}
	return m_file.m_century_table[m_file.m_century_bits.count(e)] + dist;
}

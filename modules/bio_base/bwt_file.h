#pragma once

#include "modules/bio_base/dna_base.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/io/bitcount.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/version.h"

struct bwt_header
{
	TRANSFER_OBJECT {
			VERSION(0);
			FIELD(version);
			FIELD(ca_table);
	};
	product_version version;
	std::vector<uint64_t> ca_table;
};

class bwt_range;

class bwt_file
{
	friend class bwt_range;
public:
	static constexpr uint64_t k_magic = 0x57426C6172697053; // SpiralBW
	bwt_file(const std::string& path);
	bwt_range bwt();
private:
	mmap_buffer m_buf;
	uint64_t m_offset;
	bwt_header m_header;
	const char* m_base;
	uint64_t m_entries;
	uint64_t m_bits_size;
	bitcount m_base_bits[4];
	bitcount m_century_bits;
	const uint32_t* m_century_table;

	uint64_t validate_offset(const char* buf);
};


class bwt_range
{
	friend class bwt_file;
public:
	bwt_range(const bwt_file& file, size_t begin, size_t end)
		: m_file(file), m_begin(begin), m_end(end)
	{}

	bwt_range find(const dna_sequence& s) const;
	bwt_range find(const std::string& query) const;

	bwt_range push_front(const dna_base& b) const;
	bool valid() const { return m_begin < m_end; }
	size_t begin() const { return m_begin; }
	size_t end() const { return m_end; }
	size_t matches() const { return valid() ? m_end - m_begin : 0; }
	uint32_t get_match(size_t which) const;

private:
	const bwt_file& m_file;
	size_t m_begin;
	size_t m_end;
};

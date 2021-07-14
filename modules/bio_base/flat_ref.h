
#pragma once

#include <memory>
#include <unordered_set>
#include <mutex>
#include <boost/functional/hash.hpp>

#include "modules/io/io.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/hash_io.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/mem_io.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/spec_headers.h"

class flat_ref_builder;

class flat_ref
{
	struct fixup_t;
	friend class flat_ref_builder;
	friend class spec_reader;
	friend class spec_writer;
	friend class flat_ref_tester; // With friends like this, who needs enemies?

public:
	struct index_t;
	struct extent_t;
	using extent_iter_t = std::vector<extent_t>::const_iterator;

public:
	flat_ref(const std::string& flat_file_path);
	flat_ref(std::unique_ptr<index_t> ref_index, mem_io&& raw_dna_buffer);
	~flat_ref();

	int64_t flatten(const std::string& scaffold_name, size_t pos, bool is_unknown_scaffold_fatal = true) const;
	//std::pair<std::string, size_t> expand(size_t flat_pos) const;
	dna_const_iterator get_dna(size_t flat_pos) const;
	void make_fasta(writable& fasta_writable, unsigned line_length = 80) const;

	std::vector<spec_header::scaffold_t> copy_scaffold_table() const;
	const index_t& get_index() const { return *m_index; };

	static const std::string k_magic_header;
	static std::string trim_scaffold_name(const std::string& scaffold_name);
private:
	mmap_buffer m_mmap;
	mem_io		m_mem_io;
	const char* m_dna_buf;
	size_t		m_dna_buf_size;
	std::unique_ptr<index_t> m_index;
	std::unordered_map<std::string, size_t> m_scaffold_by_name;
	mutable std::mutex m_mutex;
	mutable std::unordered_set<std::string> m_missing_contigs;

	void print_extent_fasta(writable& fasta_writable, extent_iter_t extent_iter, unsigned line_length) const;
	const char* get_raw_dna() const { return m_dna_buf; }
	size_t get_raw_dna_size() const { return m_dna_buf_size; }
	void build_scaffold_by_name();
	void build_ref_from_fasta(const std::string& flat_file_path);
	void build_ref_from_spec(const std::string& flat_file_path);
};

class flat_ref_builder
{
public:
	flat_ref_builder(readable& fasta, writable& flat);
	flat_ref_builder(readable& fasta, flat_ref& a_flat_ref);
	~flat_ref_builder();
	void run();
	void build_into_flat_ref(flat_ref& the_reference);

private:
	void start_scaffold(const std::string& name);
	void finish_scaffold();
	void add_base(char b);
	void write_base(char b);
	void finalize();
	void finish_extent();
	void build_dna_buffer();

	readable& m_fasta;
	writable& m_flat;
	std::unique_ptr<flat_ref::index_t> m_index;
	std::vector<size_t> m_n_locations;
	md5_hash_writer m_hasher;
	std::string m_scaffold_name;
	mem_io m_dna_buffer_writable;
	size_t m_cur_scaffold;
	size_t m_scaffold_offset;
	size_t m_flat_offset;
	size_t m_extent_size;
	uint8_t m_cur_byte;
	uint8_t m_cur_count;
};


// Note any deviations from the ACGTN alphabets here.
struct flat_ref::fixup_t {
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(extent_index);
		FIELD(base_index);
		FIELD(original_base);
	}

	size_t extent_index = 0;
	size_t base_index = 0;
	char original_base = 0;

	fixup_t() {}
	fixup_t(size_t the_extent_index, size_t the_base_index, char the_original_base)
		: extent_index(the_extent_index)
		, base_index(the_base_index)
		, original_base(the_original_base)
		{}
};


struct flat_ref::index_t {
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(scaffolds);
		FIELD(extents);
		FIELD(fixups);
	}

	index_t() = default;
	index_t(const index_t&) = delete;  // Move only semantics
	index_t(index_t&&) = default;

	index_t& operator=(const index_t&) = delete;  // Move only semantics
	index_t& operator=(index_t&&) = default;

	~index_t() = default;

	std::vector<spec_header::scaffold_t> scaffolds;
	std::vector<flat_ref::extent_t> extents;

	using fixup_key_t = std::pair<decltype(flat_ref::fixup_t::extent_index), decltype(flat_ref::fixup_t::base_index)>;
	using fixup_mapped_t = decltype(flat_ref::fixup_t::original_base);
	std::unordered_map<fixup_key_t, fixup_mapped_t, boost::hash<fixup_key_t>> fixups;

	static flat_ref::index_t::fixup_key_t make_fixup_key(size_t extent_index, size_t base_index)
	{
		return std::make_pair(extent_index, base_index);
	}
};

struct flat_ref::extent_t {
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(scaffold_name);
		FIELD(offset);
		FIELD(size);
		FIELD(flat);
	}

	size_t scaffold_name;   // Which logical scaffold is this part of
	size_t offset; // What is the offset of the extent in the scaffold
	size_t size;   // What is the size of the extent
	size_t flat;   // What is the offset of the extent in the flat genome

	bool operator<(const flat_ref::extent_t& rhs) const {
		if (scaffold_name != rhs.scaffold_name) {
			return scaffold_name < rhs.scaffold_name;
		}
		return offset < rhs.offset;
	}
};

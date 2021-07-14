#include <unordered_map>
#include <type_traits>

#include "base/base.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/make_unique.h"
#include "modules/io/log.h"
#include "modules/bio_format/fasta.h"
#include "modules/bio_base/flat_ref.h"
#include <endian.h>
#include <sys/mman.h>
#include <iostream>
#include <boost/algorithm/string.hpp>

const std::string flat_ref::k_magic_header = "SPR000";

flat_ref::flat_ref(const std::string& flat_file_path)
	: m_mmap(flat_file_path, mmap_buffer::mode::read_populate)
    , m_mem_io("",track_alloc("flat_ref"))
	, m_index(make_unique<index_t>())
{
	if (*m_mmap.buffer() == '>') {
		build_ref_from_fasta(flat_file_path);
	} else {
		build_ref_from_spec(flat_file_path);
	}
}

flat_ref::flat_ref(std::unique_ptr<index_t> ref_index, mem_io&& raw_dna_buffer)
	: m_mem_io(raw_dna_buffer)
	, m_dna_buf(m_mem_io.buffer())
	, m_dna_buf_size(m_mem_io.size())
	, m_index(std::move(ref_index))
{
	CHECK(m_mmap.size() == 0);
	build_scaffold_by_name();
}

flat_ref::~flat_ref()
{}

int64_t flat_ref::flatten(const std::string& scaffold_name, size_t pos, bool is_unknown_scaffold_fatal) const
{
	extent_t key;
	auto scaffold_it = m_scaffold_by_name.find(scaffold_name);
	if (scaffold_it == m_scaffold_by_name.end()) {
		if (is_unknown_scaffold_fatal) {
			throw io_exception(boost::format("Contig \"%1%\" was not found in the reference. "
				"Please use the correct reference or add \"--no-match-reference\" to override")
				% scaffold_name);
		} else {
			m_mutex.lock();
			if (!m_missing_contigs.count(scaffold_name)) {
				std::cerr << "Warning: contig '" << scaffold_name.c_str() << "' was not found in reference.\n";
				m_missing_contigs.insert(scaffold_name);
			}
			m_mutex.unlock();
			return -1;
		}
	}
	key.scaffold_name = scaffold_it->second;
	key.offset = pos;
	auto it = upper_bound(m_index->extents.begin(), m_index->extents.end(), key);
	if (it == m_index->extents.begin()) {
		return -1;
	}
	it--;
	if (it->scaffold_name != key.scaffold_name || it->offset > pos || pos >= it->offset + it->size) {
		return -1;
	}
	return it->flat + (pos - it->offset);
}

dna_const_iterator flat_ref::get_dna(size_t pos) const
{
	//TODO: Add a dna_const_iterator const ptr constructor so the const_cast can go away.
	return dna_const_iterator(const_cast<unsigned char*>(reinterpret_cast<const unsigned char*>(m_dna_buf)), pos, false);
}

void flat_ref::make_fasta(writable& fasta_writable, unsigned line_length) const
{
	CHECK(m_index);

	for (unsigned i = 0; i < m_index->scaffolds.size(); i++) {
		fasta_writable.print(">%s\n", m_index->scaffolds[i].name.c_str());

		extent_t current_scaffold_extent;
		current_scaffold_extent.scaffold_name = i;
		const auto extent_range = std::equal_range(
			m_index->extents.cbegin(),
			m_index->extents.cend(),
			current_scaffold_extent,
			[](const flat_ref::extent_t& lhs, const flat_ref::extent_t& rhs) { return lhs.scaffold_name < rhs.scaffold_name; }
		);
		for (auto extent_iter = extent_range.first; extent_iter != extent_range.second; ++extent_iter) {
			print_extent_fasta(fasta_writable, extent_iter, line_length);
		}
	}
}

void flat_ref::print_extent_fasta(writable& fasta_writable, extent_iter_t extent_iter, unsigned line_length) const
{
	unsigned long dna_base_index = 0;
	while (dna_base_index < extent_iter->size) {
		unsigned current_line_length = std::min(static_cast<unsigned long>(line_length), extent_iter->size - dna_base_index);
		dna_const_iterator dna_iter = get_dna(extent_iter->flat + dna_base_index);
		while (current_line_length--) {
			char the_base = static_cast<char>(*dna_iter++);
			auto fixup_iter = m_index->fixups
				.find(flat_ref::index_t::make_fixup_key(extent_iter - m_index->extents.cbegin(), dna_base_index));
			if (fixup_iter != m_index->fixups.end()) {
				the_base = fixup_iter->second;
			}
			fasta_writable.write(&the_base, 1);
			dna_base_index++;
		}
		fasta_writable.write("\n", 1);
	}
}

std::vector<spec_header::scaffold_t> flat_ref::copy_scaffold_table() const
{
	CHECK(m_index);
	return m_index->scaffolds;
}


void flat_ref::build_scaffold_by_name()
{
	CHECK(m_index);
	for(size_t i = 0; i < m_index->scaffolds.size(); i++) {
		std::string scaf = trim_scaffold_name(m_index->scaffolds[i].name);
		if(m_scaffold_by_name.count(scaf) > 0) {
			throw io_exception(printstring("Invalid reference: sequence identifier '%s' is not unique.",scaf.c_str()));
		}
		m_scaffold_by_name[scaf] = i;
	}
}

void flat_ref::build_ref_from_fasta(const std::string& flat_file_path)
{
	m_mmap.close();
	file_reader fasta_reader(flat_file_path);
	flat_ref_builder the_ref_builder(fasta_reader, *this);
	build_scaffold_by_name();
}

void flat_ref::build_ref_from_spec(const std::string& flat_file_path)
{
	if (m_mmap.size() < 8) {
		throw io_exception("Spec reference is way too small (< 8 bytes)");
	}

	if (std::memcmp(k_magic_header.data(), m_mmap.buffer(), k_magic_header.size()) != 0) {
		throw io_exception("The file passed as the reference is neither a Spec Reference nor a FASTA file");
	}

	const char* last_4 = m_mmap.buffer() + m_mmap.size() - 4;
	uint32_t msgpack_size = be32toh(*reinterpret_cast<const uint32_t*>(last_4));
	if (m_mmap.size() < 4 + msgpack_size) {
		throw io_exception("flat_ref reference header size too big");
	}
	const char* idx_start = m_mmap.buffer() + m_mmap.size() - (4 + msgpack_size);
	msgpack_deserialize(*m_index, std::string(idx_start, msgpack_size));
	build_scaffold_by_name();
	m_dna_buf = m_mmap.buffer() + k_magic_header.size();
	m_dna_buf_size = idx_start - (m_mmap.buffer() + k_magic_header.size());
	CHECK(m_mem_io.size() == 0);

	// advise that we'll be using the whole mmap
	madvise((void *)m_dna_buf, m_mmap.size(), MADV_WILLNEED);
}

std::string flat_ref::trim_scaffold_name(const std::string& scaffold_name)
{
	const std::string white_space_chars = " \t";
	auto first_non_white_space = scaffold_name.find_first_not_of(white_space_chars);
	auto second_white_space = scaffold_name.find_first_of(white_space_chars, first_non_white_space);
	CHECK_LT(first_non_white_space, second_white_space);

	if (first_non_white_space == std::string::npos) {
		throw io_exception(boost::format("Contig \"%1%\" contains only white space")
			% scaffold_name);
	}

	CHECK(scaffold_name.substr(first_non_white_space, second_white_space).find_first_not_of(white_space_chars) != std::string::npos);
	return scaffold_name.substr(first_non_white_space, second_white_space);
}

flat_ref_builder::flat_ref_builder(readable& fasta, writable& flat)
	: m_fasta(fasta)
	, m_flat(flat)
	, m_index(make_unique<flat_ref::index_t>())
    , m_dna_buffer_writable("",track_alloc("flat_ref:dna_buf"))
	, m_cur_scaffold(0)
	, m_scaffold_offset(0)
	, m_flat_offset(0)
	, m_extent_size(0)
	, m_cur_byte(0)
	, m_cur_count(0)
{}

flat_ref_builder::flat_ref_builder(readable& fasta, flat_ref& the_reference)
	: m_fasta(fasta)
	, m_flat(m_dna_buffer_writable)
	, m_index(make_unique<flat_ref::index_t>())
    , m_dna_buffer_writable("",track_alloc("flat_ref:dna_buf"))
	, m_cur_scaffold(0)
	, m_scaffold_offset(0)
	, m_flat_offset(0)
	, m_extent_size(0)
	, m_cur_byte(0)
	, m_cur_count(0)
{
	build_dna_buffer();
	the_reference.m_mem_io = m_dna_buffer_writable;
	the_reference.m_index = std::move(m_index);
	the_reference.m_dna_buf = the_reference.m_mem_io.buffer();
	the_reference.m_dna_buf_size = the_reference.m_mem_io.size();
}

flat_ref_builder::~flat_ref_builder() {}

void flat_ref_builder::run()
{
	m_flat.write(flat_ref::k_magic_header.data(), flat_ref::k_magic_header.size());

	build_dna_buffer();

	finalize();
}

void flat_ref_builder::build_dna_buffer()
{
	std::string line;
	bool scaffold_started = false;
	while(m_fasta.readline(line, k_maxline))
	{
		if (line[0] == '>')
		{
			if (scaffold_started) {
				finish_scaffold();
			}
			start_scaffold(boost::algorithm::trim_copy(line.substr(1)));
			scaffold_started = true;
		}
		else
		{
			for(size_t i = 0; i < line.size(); i++)
			{
				add_base(line[i]);
			}
		}
	}
	if (scaffold_started) {
		finish_scaffold();
	}

	while(m_cur_count > 0) {
		write_base('A');
	}
}

void flat_ref_builder::start_scaffold(const std::string& name)
{
	m_scaffold_name = name;
}

void flat_ref_builder::finish_scaffold()
{
	if (m_extent_size > 0) {
		finish_extent();
	}
	spec_header::scaffold_t scaf;
	scaf.name = m_scaffold_name;
	m_hasher.finish();
	scaf.md5 = m_hasher.hex();
	m_hasher.reset();
	scaf.size = m_scaffold_offset;
	m_index->scaffolds.push_back(scaf);
	m_scaffold_offset = 0;
	m_cur_scaffold++;
}

void flat_ref_builder::finish_extent()
{
	flat_ref::extent_t extent;
	extent.scaffold_name = m_cur_scaffold;
	extent.offset = m_scaffold_offset;
	extent.size = m_extent_size;
	extent.flat = m_flat_offset;
	m_index->extents.push_back(extent);

	m_scaffold_offset += m_extent_size;
	m_flat_offset += m_extent_size;
	m_scaffold_offset += m_n_locations.size();
	m_n_locations.clear();
	m_extent_size = 0;
}

void flat_ref_builder::write_base(char c)
{
	m_extent_size++;
	m_cur_byte <<= 2;
	m_cur_byte |= (int) dna_base(c);
	m_cur_count++;
	if (m_cur_count == 4) {
		m_flat.write((const char*) &m_cur_byte, 1);
		m_cur_count = 0;
		m_cur_byte = 0;
	}
};

void flat_ref_builder::add_base(char c)
{
	m_hasher.write(&c, 1);
	const char translate[]   = "AICI..GI..I.IN...IITTII.I.";
	//                         "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	const char iupac_ambig[] = ".C.A...A..G.A....AG..AA.C.";
	if (!isalpha(c)) {
		throw io_exception(printstring("Non-alpha fasta: ascii value of %d", c));
	}
	char translated = translate[toupper(c) - 'A'];
	if (translated == '.') {
		throw io_exception((boost::format("Invalid base in fasta: ascii value of %d") % c).str());
	}
	if (translated == 'N') {
		m_n_locations.push_back(m_extent_size + m_n_locations.size());
		return;
	}
	if (translated == 'I') {
		m_index->fixups.emplace(std::make_pair(flat_ref::index_t::make_fixup_key(m_index->extents.size(), m_extent_size), c));
		translated = iupac_ambig[toupper(c) - 'A'];
		CHECK_NE(translated, '.');
	}

	if (m_n_locations.size() > 10) {
		finish_extent();
	}
	if (m_n_locations.size() > 0) {
		while (m_n_locations.size()) {
			write_base('A');
			m_index->fixups.emplace(std::make_pair(flat_ref::index_t::make_fixup_key(m_index->extents.size(), m_n_locations.back()), 'N'));
			m_n_locations.pop_back();
		}
		CHECK(m_n_locations.empty());
	}
	write_base(translated);
}

void flat_ref_builder::finalize()
{
	while(m_cur_count > 0) {
		write_base('A');
	}
	std::string s = msgpack_serialize(*m_index);
	m_flat.write(s.data(), s.size());
	uint32_t no_size = htobe32(uint32_t(s.size()));
	m_flat.write(reinterpret_cast<char*>(&no_size), 4);
	m_flat.close();
}

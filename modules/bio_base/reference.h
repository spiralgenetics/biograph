#pragma once

#include "modules/bio_base/bwt_file.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/flat_ref.h"
#include "modules/bio_base/reference_assembly.h"
#include "modules/bio_base/seq_position.h"
#include "modules/io/mmap_buffer.h"

class reference
{
public:
	reference(const std::string& ref_name);
	reference(const std::string& ref_name, const std::string& reference_assembly_json_parent_dir);

	const reference_assembly& get_assembly() const { return m_reference_assembly; }
	const std::string& ref_name() const { return m_ref_name; }
    
	seq_position get_seq_position(size_t pos) const { return m_reference_assembly.get_seq_position(pos); }
	seq_position get_seq_position(dna_const_iterator it) const;
	size_t flatten(const seq_position& pos) const { return m_reference_assembly.flatten(pos); }
	size_t flatten(std::string scaffold_name, size_t pos) const { return m_reference_assembly.flatten(scaffold_name, pos); }

	std::pair<size_t, size_t> flatten_range(
		const std::string& contig_name
		, unsigned long start
		, unsigned long end
		, bool use_exact_loci = true
	) const
	{ return m_reference_assembly.flatten_range(contig_name, start, end, use_exact_loci); }

	dna_const_iterator get_dna(size_t pos) const;
	dna_slice get_supercontig(size_t pos) const;
	size_t size() const { return m_reference_assembly.size(); }
	bwt_range get_bwt() const;
	const flat_ref& get_flat_ref() const;

	std::string path() const { return m_path; }
	std::string fasta_path() const;

private:
	std::string m_ref_name;
	std::string m_path;
	reference_assembly m_reference_assembly;
	mutable mmap_buffer m_dna_buf;

	// mutable refs for lazy loading
	mutable std::unique_ptr<flat_ref> m_flat_ref;
	mutable std::unique_ptr<bwt_file> m_bwt;

	// these must be member variables. Otherwise unittest will fail to load
	// the resources on subsequent tests.
	mutable std::once_flag init_ref_once_flag;
	mutable std::once_flag init_bwt_once_flag;
};

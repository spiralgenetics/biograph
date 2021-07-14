#include <mutex>

#include <boost/algorithm/string.hpp>

#include "modules/bio_base/bwt_file.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/karyotype_compat.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/defaults.h"
#include "modules/io/make_unique.h"

reference::reference(const std::string& ref_name)
	: reference(ref_name, CONF_S(reference_path))
{}

reference::reference(const std::string& ref_name, const std::string& reference_assembly_json_dir)
	: m_ref_name(ref_name)
	, m_path(reference_assembly_json_dir + "/" + ref_name + "/")
{
	using namespace boost;
	
	// Compress trailing //s to one
	m_path = trim_right_copy_if(m_path, is_any_of("/")) + "/";

	// import supercontigs to m_reference_assembly
	kt_compat m_karyotype;
	json_deserialize(m_karyotype, slurp_file(m_path + defaults.karyotype));

	m_reference_assembly.scaffolds = m_karyotype.chromosomes;
	m_reference_assembly.scaffold_order = m_karyotype.chr_order;

	for(const auto& ktsc : m_karyotype.supercontigs) {
		supercontig sc(ktsc.chr, ktsc.offset, ktsc.len);
		m_reference_assembly.supercontigs.insert(sc);
	}

	m_reference_assembly.generate_tables();
}

dna_const_iterator reference::get_dna(size_t pos) const
{
	if(!m_flat_ref) {
		std::call_once(init_ref_once_flag, [this]()
		{
			m_flat_ref = make_unique<flat_ref>(m_path + defaults.reference_ref);
		});
	}
	return m_flat_ref->get_dna(pos);
}

bwt_range reference::get_bwt() const
{
	if(!m_bwt) {
		std::call_once(init_bwt_once_flag, [this]()
		{
			m_bwt = make_unique<bwt_file>(m_path + defaults.reference_bwt);
		});
	}
	return m_bwt->bwt();
}

const flat_ref& reference::get_flat_ref() const
{
  if (!m_flat_ref) {
    std::call_once(init_ref_once_flag, [this]() {
			m_flat_ref = make_unique<flat_ref>(m_path + defaults.reference_ref);
		});
  }
  return *m_flat_ref;
}

dna_slice reference::get_supercontig(size_t pos) const
{
	const supercontig& sc = m_reference_assembly.get_supercontig(pos);
	return dna_slice(get_dna(sc.tot_offset), get_dna(sc.tot_offset + sc.len));
}

seq_position reference::get_seq_position(dna_const_iterator it) const
{
	return get_seq_position((size_t) it.get_original_offset());
}

std::string reference::fasta_path() const
{
	return m_path + defaults.reference_fasta;
}


#pragma once

#include <string>
#include <boost/filesystem.hpp>
// Default constants go here.
class Defaults {
	public:
		const std::string original_fasta = "source.fasta";
		const std::string reference_fasta = "reference.fasta";
		const std::string reference_ref = "reference.ref";
		const std::string reference_seqset = "reference.seqset";
		const std::string reference_bwt = "reference.bwt";
		const std::string alu_fasta = "alu.fasta";
		const std::string karyotype = "karyotype.json";

		//See if each of the ref files are there like they should be
		bool check_refdir(std::string& refdir) {

			namespace fs = boost::filesystem;

			fs::path m_path(refdir);

			if (!fs::is_directory(refdir)) {
				std::cerr << "--ref argument  '" << refdir << "' is not a directory. ";
				std::cerr << "Note that FASTA files must be prepared using 'biograph reference'." << std::endl;
				return false;
			}

			if (!fs::exists(m_path.append(original_fasta))){
				std::cerr << original_fasta << " missing from reference directory " << refdir << std::endl;
				return false;
			}
			m_path = m_path.parent_path();

			if (!fs::exists(m_path.append(reference_fasta))){
				std::cerr << m_path.append(reference_fasta) << std::endl;
				std::cerr << reference_fasta << " missing from reference directory " << refdir << std::endl;
				return false;
			}
			m_path = m_path.parent_path();

			if (!fs::exists(m_path.append(reference_ref))){
				std::cerr << reference_ref << " missing from reference directory " << refdir << std::endl;
				return false;
			}
			m_path = m_path.parent_path();
			/* Apparently this isn't made, yet or used anywhere...
			if (!fs::exists(m_path.append(reference_seqset))){
				std::cerr << reference_seqset << " missing from directory " << refdir << std::endl;
				return false;
			}
			m_path = m_path.parent_path();
			*/

			if (!fs::exists(m_path.append(reference_bwt))){
				std::cerr << reference_seqset << " missing from reference directory " << refdir << std::endl;
				return false;
			}
			m_path = m_path.parent_path();
			/* Apparently optional
			if (!fs::exists(m_path.append(alu_fasta))){
				std::cerr << alu_fasta << " missing from directory " << refdir << std::endl;
				return false;
			}
			m_path = m_path.parent_path();
			*/

			if (!fs::exists(m_path.append(karyotype))){
				std::cerr << karyotype << " missing from reference directory " << refdir << std::endl;
				return false;
			}
			m_path = m_path.parent_path();

			return true;
		}
};

static Defaults defaults;

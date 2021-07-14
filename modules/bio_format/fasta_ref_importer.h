#pragma once

#include "modules/bio_base/reference_assembly.h"
#include "modules/mapred/path.h"
#include "modules/io/file_io.h"
#include "modules/io/progress_tracker_types.h"
#include "modules/io/progress_tracker.h"
#include "modules/bio_format/importer.h"

struct supercontig_datum
{
	unsigned long m_index;
	std::string m_name;
	std::string m_sequence;
	unsigned long m_start;

	supercontig_datum(unsigned long index, const std::string& name, const std::string& sequence, unsigned long start)
		: m_index(index), m_name(name), m_sequence(sequence), m_start(start) {}

	bool operator<(const supercontig_datum& right) const
	{
		if (m_index < right.m_index) return true;
		else if (m_index == right.m_index) return m_start < right.m_start;
		else return false;
	}
};

class fasta_ref_importer
{
public:
  fasta_ref_importer(const path& out_dir, readable& in, const std::vector<std::string>& scaffold_order, size_t min_n_run, progress_t& update);
	void run();
	const reference_assembly& get_assembly() { return m_reference_assembly; }

private:
	void add_base(char c);
	void finish_supercontig();
	void finish_scaffold();
	void store_supercontig();
	void write_supercontigs();

	reference_assembly m_reference_assembly;
	path m_out_dir;
	path m_seq_subdir;
  size_t m_min_n_run;
	readable& m_fasta_in;
	std::string m_scaffold_name;  // Current scaffold name
	size_t m_position;  // Current position
	size_t m_start;  // Start of sequence
	size_t m_end;  // End of sequence
	bool m_had_contig;  // Used to remove 'empty' scaffolds
	std::string m_sequence;  // Actual in progress sequence as string
	std::unique_ptr<writable> m_fasta_out;
	std::vector<supercontig_datum> m_supercontig_data;
	progress_tracker m_tracker;
	size_t m_total_bytes_read;
};


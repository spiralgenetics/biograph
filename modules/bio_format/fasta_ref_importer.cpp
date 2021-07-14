#include "modules/bio_format/fasta_ref_importer.h"
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <regex>
#include "modules/bio_base/karyotype_compat.h"
#include "modules/bio_format/fasta.h"
#include "modules/io/defaults.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"

Defaults spiral_defaults;

fasta_ref_importer::fasta_ref_importer(const path& out_dir, readable& in,
                                       const std::vector<std::string>& scaffold_order,
                                       size_t min_n_run,
                                       progress_t& update)
    : m_out_dir(out_dir),
      m_min_n_run(min_n_run),
      m_fasta_in(in),
      m_fasta_out(out_dir.append(spiral_defaults.reference_fasta).write()),
      m_tracker(update),
      m_total_bytes_read(0) {
  m_reference_assembly.scaffold_order = scaffold_order;
  m_supercontig_data.reserve(scaffold_order.size());
}

void fasta_ref_importer::run() {
  const char* prefix = "fasta_ref_importer::run>";
  SPLOG("%s begin", prefix);

  std::string line;
  boost::regex name_extractor("^>(\\S+)\\s*.*$");
  while (m_fasta_in.readline(line, k_maxline)) {
    if (line[0] == '>') {
      finish_scaffold();
      boost::smatch base_match;
      if (boost::regex_match(line, base_match, name_extractor)) {
        m_scaffold_name = base_match[1];
      } else {
        throw io_exception(boost::str(boost::format("Illegal fasta entry name '%s'") % line.c_str()));
      }
    } else {
      for (size_t i = 0; i < line.size(); i++) {
        add_base(line[i]);
        m_total_bytes_read++;
        m_tracker.update(m_total_bytes_read, m_total_bytes_read);
      }
    }
  }
  SPLOG("%s done reading the fasta file", prefix);

  finish_scaffold();

  // export m_reference_assembly to supercontigs

  kt_compat karyotype;

  karyotype.chromosomes = m_reference_assembly.scaffolds;
  karyotype.chr_order = m_reference_assembly.scaffold_order;

  for (const auto& sc : m_reference_assembly.supercontigs) {
    kt_supercontig ktsc(sc.scaffold_name, sc.offset, sc.len);
    karyotype.supercontigs.insert(ktsc);
  }

  m_out_dir.append(spiral_defaults.karyotype).put(json_serialize(karyotype));
  write_supercontigs();

  SPLOG("%s end", prefix);
}

void fasta_ref_importer::add_base(char c) {
  const char translate[] = "ANCN..GN..N.NN...NNTTNNNN.";
  if (!isalpha(c)) throw io_exception(printstring("Non-alpha fasta: ascii value of %d", c));
  c = translate[toupper(c) - 'A'];
  if (c == '.')
    throw io_exception((boost::format("Invalid base in fasta: ascii value of %d") % c).str());

  if (c != 'N') {
    if (m_end != m_position) {
      // There was a break (some N's).  Should we start a new sequence or
      // just add the N's in?  We keep them connected in the number of N's is < 50
      if (m_position - m_end < m_min_n_run) {
        while (m_end != m_position) {
          m_sequence.push_back(c);
          m_end++;
        }
      } else
        finish_supercontig();
    }
    m_sequence.push_back(c);
    m_end++;
  }
  m_position++;
}

void fasta_ref_importer::finish_supercontig() {
  if (m_sequence.size() == 0) {
    m_start = m_end = m_position;
    return;
  }
  m_had_contig = true;

  store_supercontig();

  // Add to reference_assembly
  m_reference_assembly.supercontigs.insert(
      supercontig(m_scaffold_name, m_start, m_sequence.size()));

  // Reset state
  m_start = m_end = m_position;
  m_sequence.clear();
}

void fasta_ref_importer::store_supercontig() {
  auto seq_name_iter = std::find(m_reference_assembly.scaffold_order.begin(),
                                 m_reference_assembly.scaffold_order.end(), m_scaffold_name);
  auto scaffold_index = std::distance(m_reference_assembly.scaffold_order.begin(), seq_name_iter);
  supercontig_datum a_datum(scaffold_index, m_scaffold_name, m_sequence, m_start);
  m_supercontig_data.push_back(a_datum);
}

void fasta_ref_importer::finish_scaffold() {
  if (m_scaffold_name != "") {
    finish_supercontig();
    if (m_had_contig) {
      int index = -1;
      for (int i = 0; i < (int)m_reference_assembly.scaffold_order.size(); i++) {
        if (m_reference_assembly.scaffold_order[i] == m_scaffold_name) {
          index = i;
          break;
        }
      }
      if (index == -1) {
        index = m_reference_assembly.scaffold_order.size();
        m_reference_assembly.scaffold_order.push_back(m_scaffold_name);
      }

      m_reference_assembly.scaffolds.insert(scaffold(m_scaffold_name, m_position, index));
    }
  }
  m_position = m_start = m_end = 0;
  m_had_contig = false;
}

void fasta_ref_importer::write_supercontigs() {
  const char* prefix = "fasta_ref_importer::write_supercontigs>";
  SPLOG("%s begin", prefix);

  std::sort(m_supercontig_data.begin(), m_supercontig_data.end());

  // Write to fasta output file
  for (unsigned int i = 0; i < m_supercontig_data.size(); i++) {
    m_fasta_out->print(">%s:%ld", m_supercontig_data[i].m_name.c_str(),
                       m_supercontig_data[i].m_start);
    for (unsigned int j = 0; j < m_supercontig_data[i].m_sequence.size(); j++) {
      if (j % 80 == 0) {
        m_fasta_out->print("\n");
      }
      m_fasta_out->print("%c", m_supercontig_data[i].m_sequence[j]);
      m_total_bytes_read++;
      m_tracker.update(m_total_bytes_read, m_total_bytes_read);
    }
    m_fasta_out->print("\n");
  }
  SPLOG("%s end", prefix);
}

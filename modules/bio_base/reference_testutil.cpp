#include "modules/bio_base/reference_testutil.h"
#include "modules/test/build_ref.h"
#include "modules/io/make_unique.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include <iostream>
#include <gtest/gtest.h>

std::unique_ptr<flat_ref> create_flat_ref(std::vector<dna_sequence> seqs) {
  dna_sequence all_seqs;
  std::unique_ptr<flat_ref::index_t> index(new flat_ref::index_t);

  size_t scaffold_idx = 0;
  for (const dna_sequence& seq : seqs) {
    flat_ref::extent_t extent = {};
    extent.scaffold_name = scaffold_idx;
    extent.offset = 0;
    extent.size = seq.size();
    // Apparently this "flat" is 1-indexed instead of 0-indexed?
    // TODO(nils): Investigate this further and document more clearly.
    extent.flat = all_seqs.size() + 1;
    index->extents.push_back(extent);

    spec_header::scaffold_t scaffold = {};
    scaffold.name = std::to_string(scaffold_idx);
    scaffold.size = extent.size;
    index->scaffolds.push_back(scaffold);

    all_seqs += seq;
    ++scaffold_idx;
  }

  mem_io raw_dna_buffer(all_seqs.as_packed(), track_alloc("reference_testutil:flat_ref_dna"));
  return make_unique<flat_ref>(std::move(index), std::move(raw_dna_buffer));
}

std::unique_ptr<reference> create_reference(
    const std::vector<dna_sequence>& seqs) {
  std::vector<std::string> seq_strings;
  for (const auto& seq : seqs) {
    seq_strings.push_back(seq.as_string());
  }
  return create_reference_str(seq_strings);
}

std::unique_ptr<reference> create_reference_str(const std::vector<std::string>& seqs) {
  static size_t reference_num = 0;
  size_t this_reference_num = reference_num++;

  boost::filesystem::create_directories(CONF_S(temp_root));
  boost::filesystem::create_directories(CONF_S(reference_path));

  std::string fasta_path =
      CONF_S(temp_root) +
      "/ref" + std::to_string(this_reference_num) + ".fasta";
  {
    file_writer f(fasta_path);

    size_t scaffold_num = 0;
    for (const std::string& seq : seqs) {
      size_t sc = scaffold_num++;
      std::string record =
          ">" + std::to_string(sc) + "\n" + seq + "\n";
      f.write(record.data(), record.size());
    }

    // For test data, also insert 4 scaffolds containing the
    // individual bases so we don't have to make sure we have at least
    // 1 of each base:
    for (dna_base b : dna_bases()) {
      size_t sc = scaffold_num++;
      std::string record =
          ">" + std::to_string(sc) + "\n" + char(b) + "\n";
      f.write(record.data(), record.size());
    }
  }

  std::string ref_name = "ref" + std::to_string(this_reference_num);
  perform_build_ref(ref_name, fasta_path);
  return make_unique<reference>(ref_name);
}

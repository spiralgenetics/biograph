#pragma once
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/reference.h"
#include "python/biograph/seqset.h"
#include <pybind11/pybind11.h>

// A python reference is a lightweight reference locus, which translates to C++
// as a reference object on the heap accessed by a shared_ptr and a pair of
// seq_positions into that reference.  The seq positions in python are string-
// unsigned long pairs where the string is the contig name rather than the index
// into the reference assembly.

typedef pybind11::list py_list;
typedef pybind11::dict py_dict;

class reference_range;

// Python wrapper for bwt manipulation
class bwt_wrapper {
 public:
  // Construct a new wrapper with the bwt_range pointing to the location of the
  // sequence.
  bwt_wrapper(const std::shared_ptr<reference> the_reference,
              const dna_sequence& seq)
      : m_reference(the_reference),
        m_bwt(the_reference->get_bwt().find(seq)),
        m_query(seq) {}

  std::string str() const;
  std::string repr() const;

  bwt_wrapper find(const dna_sequence& s);

  bool valid() const;
  size_t begin() const;
  size_t end() const;
  size_t matches() const;

  std::string ref_path() const { return m_reference->path(); }

  // Return a reference_range for the given match
  // throw (raise a python RuntimeError) if the bwt is not valid
  reference_range get_match(size_t which) const;

  dna_sequence query() const { return m_query; };

 private:
  std::shared_ptr<reference> m_reference;
  bwt_range m_bwt;
  dna_sequence m_query;
};

// reference ranges
class reference_range {
 public:
  reference_range(std::shared_ptr<reference> the_reference, size_t flat_start,
                  size_t flat_end);

  std::string str() const;
  std::string repr() const;

  std::string get_scaffold() const;
  unsigned long get_start() const;
  unsigned long get_end() const;

  size_t get_flat_start() const { return m_flat_ref_start; }
  size_t get_flat_end() const { return m_flat_ref_end; }

  dna_sequence sequence() const;
  size_t size() const;

  std::shared_ptr<reference> get_reference() const { return m_reference; }

 private:
  std::shared_ptr<reference> m_reference;
  size_t m_flat_ref_start;
  size_t m_flat_ref_end;
};

// Python wrapper for reference objects
class reference_wrapper {
 public:
  reference_wrapper(const std::string& reference_dir_path)
      : m_reference(std::make_shared<reference>("", reference_dir_path)) {}

  std::string str() const;
  std::string repr() const;

  reference_range make_range(const std::string& contig_name,
                             unsigned long start, unsigned long end,
                             bool use_exact_loci = true) const;

  size_t size() const;
  py_list scaffolds() const;
  py_list supercontigs() const;
  reference_range get_supercontig(std::string supercontig_name,
                                  size_t pos) const;

  py_list find_ranges(const std::string& supercontig_name, size_t start,
                      size_t end) const;
  py_dict scaffold_lens() const;

  std::string path() const { return m_reference->path(); }

  bwt_wrapper find_sequence(const dna_sequence& seq) const;
  bwt_wrapper find(const std::string& seq_str) const;

  std::shared_ptr<reference> get_reference() const { return m_reference; }

 private:
  std::shared_ptr<reference> m_reference;
};

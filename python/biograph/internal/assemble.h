#pragma once

#include <pybind11/pybind11.h>
#include "modules/bio_base/call_structural.h"
#include "modules/bio_base/seqset_assemble.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"

typedef pybind11::list py_list;

class variant_wrapper {
 public:
  variant_wrapper(const sv_out& var, std::shared_ptr<assembly> src,
                  std::shared_ptr<reference> ref);

  std::string str() const;
  std::string repr() const;

  bool py_lessthan(variant_wrapper other);
  int py_cmp(variant_wrapper other);

  bool is_structural() const { return m_var.is_structural; }

  std::string left_contig() const;
  std::string right_contig() const;
  unsigned long left_position() const;
  unsigned long right_position() const;

  bool left_forward() const { return !m_var.left_ref.is_rev_comp(); }
  bool right_forward() const { return !m_var.left_ref.is_rev_comp(); }

  dna_sequence variant_sequence() const;
  dna_sequence assembly_sequence() const;

  reference_range range() const;  // Only valid if in same supercontig, both fwd

  int assembly_begin() const { return m_var.seq_begin; }
  int assembly_end() const { return m_var.seq_end; }

  py_list depths() const;
  py_list assembly_depths() const;

  uint64_t assembly_id() const { return m_assembly->id; }

  uint8_t min_depth() const { return m_min_depth; }
  uint8_t max_depth() const { return m_max_depth; }
  float avg_depth() const { return m_avg_depth; }
  uint8_t min_overlap() const { return m_assembly->min_overlap; }
  variant_wrapper flip() const;

 private:
  variant_wrapper() {}
  sv_out m_var;
  std::shared_ptr<assembly> m_assembly;
  std::shared_ptr<reference> m_ref;
  uint8_t m_min_depth;
  uint8_t m_max_depth;
  float m_avg_depth;
};

py_list py_assemble(py_list src, py_list dest, uint8_t min_overlap,
                    uint32_t max_steps, bool skip_ambig,
                    const std::shared_ptr<readmap>& wrap);

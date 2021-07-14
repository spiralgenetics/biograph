#pragma once

#include <pybind11/pybind11.h>
#include "modules/bio_base/seqset_anchor.h"
#include "python/biograph/biograph.h"

typedef pybind11::list py_list;

class anchor_wrapper {
  friend py_list py_assemble(py_list src, py_list dest, uint8_t min_overlap,
                             uint32_t max_steps, bool skip_ambig,
                             const std::shared_ptr<readmap>& wrap);

 public:
  anchor_wrapper(const anchor& anchor, std::shared_ptr<seqset_file> seqset,
                 std::shared_ptr<reference> reference);
  seqset_range read() const;
  reference_range range() const;
  bool forward() const;
  uint8_t overlap() const { return m_anchor.overlap; }

 private:
  anchor m_anchor;
  std::shared_ptr<seqset_file> m_seqset;
  std::shared_ptr<reference> m_reference;
};

py_list py_anchor(const std::shared_ptr<biograph>& wrap,
                  const reference_range& ref_range, bool forward = true,
                  uint8_t min_overlap = 70, uint32_t max_anchors = 10000);

#pragma once

#include "modules/variants/assemble.h"

namespace variants {

class scaffold {
 public:
  struct extent {
    aoffset_t offset = 0;
    dna_slice sequence;
  };

  class iterator {
   public:
    iterator() = default;
    iterator(const iterator&) = default;
    iterator& operator=(const iterator&) = default;

    dna_base operator*() const { return *m_extent_it; }
    dna_const_iterator operator->() const { return m_extent_it; }
    iterator& operator++() {
      ++m_extent_it;
      ++m_offset;
      if (m_extent_it == m_scaffold_it->sequence.end()) {
        ++m_scaffold_it;
        if (m_scaffold_it == m_scaffold->extents().end()) {
          m_offset = m_scaffold->end_pos();
        } else {
          m_offset = m_scaffold_it->offset;
          m_extent_it = m_scaffold_it->sequence.begin();
        }
      }
      return *this;
    }

    aoffset_t offset() const { return m_offset; }
    bool first_in_extent() const { return m_extent_it == m_scaffold_it->sequence.begin(); }
    aoffset_t extent_end_offset() const {
      CHECK(m_scaffold_it != m_scaffold->extents().end());
      return m_scaffold_it->offset + aoffset_t(m_scaffold_it->sequence.size());
    }

    bool operator==(const iterator& rhs) const { return m_offset == rhs.m_offset; }
    bool operator!=(const iterator& rhs) const { return !((*this) == rhs); }

    // Description is debug info for where this skip comes from, in
    // case of problems.  TODO(nils): Remove "description" arg once we
    // save debug symbols for release builds.
    void skip_to(aoffset_t offset, const char* description);

   private:
    friend class scaffold;
    iterator(const scaffold* s, std::vector<scaffold::extent>::const_iterator scaffold_it,
             dna_const_iterator extent_it, aoffset_t offset)
        : m_scaffold(s), m_scaffold_it(scaffold_it), m_extent_it(extent_it), m_offset(offset) {}

    const scaffold* m_scaffold = nullptr;
    std::vector<scaffold::extent>::const_iterator m_scaffold_it;
    dna_const_iterator m_extent_it;
    aoffset_t m_offset = 0;
  };

  scaffold() = default;
  scaffold(dna_slice simple) { add(0, simple); }
  scaffold(const dna_sequence& simple) { add(0, simple); }
  scaffold(const std::vector<extent>& extents) : m_extents(extents) { m_end_pos = calc_end_pos(); }
  scaffold(const std::vector<extent>& extents, aoffset_t end_pos)
      : m_extents(extents), m_end_pos(end_pos) {
    CHECK_GE(m_end_pos, calc_end_pos());
  }

  void add(aoffset_t offset, dna_slice seq);
  void add(aoffset_t offset, const dna_sequence& seq) { add(offset, save_storage(seq)); }

  const std::vector<extent>& extents() const { return m_extents; }

  bool empty() const { return m_extents.empty(); }

  scaffold subscaffold(aoffset_t start, aoffset_t len) const;
  std::pair<dna_slice, dna_slice> split_extent_at(aoffset_t start) const;
  bool is_simple() const;
  dna_slice get_simple() const;

  aoffset_t end_pos() const { return m_end_pos; }
  void set_end_pos(aoffset_t new_end_pos);

  friend std::ostream& operator<<(std::ostream& os, const scaffold& s) {
    s.print_to(os);
    return os;
  }

  void print_to(std::ostream& os) const;

  std::string as_string() const { return subscaffold_str(0, end_pos()); }
  std::string subscaffold_str(aoffset_t start, aoffset_t len) const;

  unsigned shared_prefix_length(dna_slice seq) const;
  scaffold rev_comp() const;

  iterator begin() const;
  iterator end() const;

 private:
  aoffset_t calc_end_pos() const;
  dna_slice save_storage(const dna_sequence& seq);
  void reverse_in_place();

  std::vector<extent> m_extents;
  aoffset_t m_end_pos = 0;
  std::shared_ptr<std::list<dna_sequence>> m_seq_storage;
};

}  // namespace variants

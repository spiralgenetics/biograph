#pragma once

#include <boost/icl/interval_set.hpp>
#include <boost/icl/separate_interval_set.hpp>
#include <deque>
#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

// Splits assemblies into variants.  Assemblies that cover the same
// regions of reference will be padded.  Reference assemblies are also
// generated.  Regions of reference without any variant assembies are
// skipped.
//
// All coverage (depth) counts refer to the coverage between two bases.
// For a sequence ABCD, the bases would have indexes (0, 1, 2, 3).  The coverage would have index 0
// = coverage between A and B, index 1 = coverage between B and C, and index 2 = coverage between C
// and D.

class split_variants : public sorted_output_pipeline_step {
 public:
  split_variants(const assemble_options& options, pipeline_step_t output);
  virtual ~split_variants();

  void on_assembly(assembly_ptr a) override;

 private:
  struct coverage_entry {
    aoffset_t offset = 0;
    int depth = 0;
    unsigned min_overlap = 0;
    std::vector<uint32_t> read_ids;

    bool operator<(aoffset_t rhs) const { return offset < rhs; }
    friend bool operator<(aoffset_t rhs, const coverage_entry& e) { return rhs < e.offset; }
  };

  struct coverage_state {
    std::deque<coverage_entry> coverage;
    aoffset_t offset = 0;
    seqset_range cur;
  };

  struct asm_info {
    assembly_ptr a;
    coverage_state cov;
  };

  struct ref_info {
    scaffold::iterator scaffold_it;
    coverage_state cov;
    std::unordered_set<uint32_t, unsalted_hash> seen_read_ids;
  };

  void advance_ref_coverage_range(aoffset_t flush_to, aoffset_t target);
  assembly_ptr make_ref_coverage_assembly(aoffset_t start, aoffset_t limit);

  // TODO(nils): For assemblies that match reference in some places,
  // add their coverage to the ref coverage.
  void add_coverage_to_ref(const std::vector<coverage_entry>& entries, aoffset_t ref_start,
                           aoffset_t seq_start, aoffset_t seq_limit);

  void chunk_and_output();
  void output_variant_regions(const asm_info& i);

  void advance_coverage(coverage_state* cov, dna_base base);

  // Returns the min depth for the closed interval [start, limit]
  void get_min_depth(const coverage_state& cov, aoffset_t start, aoffset_t limit);

  assemble_options m_options;

  aoffset_t m_leftmost_offset = 0;
  aoffset_t m_rightmost_offset = 0;

  ref_info m_ref;
  boost::icl::interval_set<aoffset_t> m_variant_regions;
  std::vector<asm_info> m_active;
};
};

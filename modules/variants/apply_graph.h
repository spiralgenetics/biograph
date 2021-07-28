#pragma once

#include "modules/variants/assemble.h"

namespace variants {

struct graph_context {
  // Variant assembly in question
  assembly* a = nullptr;

  // Assemblies on the reference branch for this variant.
  std::vector<assembly*> refs;

  // Adjacent reference assemblies on either side
  assembly* left_ref = nullptr;
  assembly* right_ref = nullptr;

  read_coverage_t ref_coverage() const { return merge_coverage(&assembly::read_coverage); }
  read_coverage_t ref_pair_coverage() const {
    return merge_coverage(&assembly::pair_read_coverage);
  }

 private:
  read_coverage_t merge_coverage(boost::optional<read_coverage_t> assembly::*field) const;

  /* IDEA(nils): Add a list of overlapping assemblies also */
};

// Applies a function to each non-reference assembly, including
// information on surrounding and overlapping assemblies.
class apply_graph : public sorted_output_pipeline_step {
 public:
  using on_context_func_t = std::function<void(graph_context)>;

  apply_graph(const on_context_func_t& f, pipeline_step_t output);
  ~apply_graph() { flush(); }

  void on_assembly(assembly_ptr a) override;

  void flush() override;

 private:
  using aptr = explicit_shared_ptr<assembly_ptr, false /* no need to be atomic */,
                                   false /* do not allow implicit copy */,
                                   false /* do not allow implicit delete */>;

  struct result {
    aptr a;

    aptr left_ref;
    aptr right_ref;

    // TODO(nils): Is there a way to use aptr with std::vector without
    // it trying to invoke the copy constructor?
    std::list<aptr> refs;

    result() = default;
    result(const result&) = delete;
    result(result&&) = default;
    result& operator=(result&&) = default;
  };

  void advance_to(aoffset_t target);
  void advance_towards(aoffset_t target);
  void output_result(result r);
  void release(aptr a);

  on_context_func_t m_on_context;
  aoffset_t m_cur_offset = std::numeric_limits<aoffset_t>::min();

  aptr m_left_ref;
  aptr m_right_ref;

  absl::btree_multimap<aoffset_t /* right offset */, result> m_active;
};

}  // namespace variants

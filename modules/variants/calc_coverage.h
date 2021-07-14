#pragma once

#include "modules/variants/assemble.h"
#include "modules/variants/path_group.h"
#include "modules/variants/scaffold.h"

namespace variants {

// calc_coverage calculates interbase coverage for each incoming assembly, along with reference
// coverage.
//
// Note that interbase coverage is the coverage between each base.  An
// assembly with a sequence length of 100 will get 101 coverage
// values.  The first and last correspond to the coverage to the left
// of the first base and to the right of the last base.
class calc_coverage : public sorted_output_pipeline_step, public path_group::listener {
 public:
  calc_coverage(const assemble_options& options, pipeline_step_t output);
  calc_coverage() = delete;
  ~calc_coverage() override;

  void on_assembly(assembly_ptr a) override;

 private:
  int sum_dir_depth(int fwd, int rev) const;

  class coverage_accum {
   public:
    void initialize(int* starts_and_ends, size_t cov_size) {
      m_starts = starts_and_ends;
      m_ends = m_starts + cov_size;
    }
    static size_t cov_elems_needed(size_t cov_size) {
      return cov_size /* starts */ + cov_size /* ends */;
    }

    // Returns true if this read contributed at all to the coverage.
    bool add(int start, int end);
    std::vector<int> coverage();
    size_t size() const { return m_ends - m_starts; }

   private:
    int* m_starts = nullptr;
    int* m_ends = nullptr;
  };
  struct other_min_coverage {
    int base = 0;
    int pair = 0;
  };
  struct shared_assembly {
    assembly_ptr a;
    calc_coverage* c;
    // List of read ids to match pairs with
    std::unordered_set<uint32_t, unsalted_hash> pair_read_ids;
    // List of minimum coverages of other non-reference assemblies in the same place.
    std::vector<other_min_coverage> other_min;
    // List of coverages of ref sections in the same place; there may
    // be multiple ref sections as assemblies start and end.
    std::vector<std::shared_ptr<std::vector<int>>> other_ref_coverage;
    ~shared_assembly();
  };
  struct cov_tracker : public path_group::distant_object {
    std::shared_ptr<shared_assembly> var_assembly;
    std::vector<std::shared_ptr<shared_assembly>> other_assemblies;
    coverage_accum fwd_coverage;
    coverage_accum rev_coverage;
    coverage_accum pair_coverage;
    calc_coverage* c = nullptr;
    std::unique_ptr<int[]> starts_and_ends;
    ~cov_tracker();
  };

  class dobj_visitor : public path_group::dobj_visitor {
   public:
    dobj_visitor(const seqset_range& r, const readmap* readmap, std::pair<uint32_t, uint32_t> reads)
        : m_r(r), m_readmap(readmap), m_reads(reads) {}
    void visit(path_group::distant_object*, int distance) override;

   private:
    seqset_range m_r;
    const readmap* m_readmap = nullptr;
    std::pair<uint32_t, uint32_t> m_reads;
  };

  void on_seqset_entry(const seqset_range& r, path_group* pg) override;
  void skip_to(aoffset_t offset);
  void advance_to(aoffset_t offset);
  void advance_towards(aoffset_t offset);
  void init_ref_pg();
  void add_ref_dobj(int ref_len);
  void flush_active_to_here();

  std::multimap<aoffset_t /* right offset */,
                std::pair<std::shared_ptr<shared_assembly>, std::unique_ptr<path_group>>>
      m_active;
  std::unique_ptr<path_group> m_ref_path_group;
  assemble_options m_options;

  aoffset_t m_cur_offset = 0;
  aoffset_t m_need_coverage_to = 0;
  scaffold::iterator m_scaffold_it;
  bool m_need_ref_dobj = false;
};

}  // namespace variants

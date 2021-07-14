// A class to provide vargraph-like coverage on assemblies.
//
// Inputs are:
//   A list of assemblies
//

#include <boost/container/small_vector.hpp>

#include "absl/container/flat_hash_set.h"
#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

class read_cov : public sorted_output_pipeline_step {
 public:
  read_cov(const assemble_options& opts, pipeline_step_t output);
  ~read_cov();

  void add_variant(aoffset_t start, aoffset_t limit, const dna_sequence& seq);

  void on_assembly(assembly_ptr a) override;

  void on_path_trim(size_t paths);

 private:
  struct result;
  using roff_index_t = uint32_t;
  struct result_offset;
  struct keep_result;
  friend std::ostream& operator<<(std::ostream& os, const read_cov::result_offset& roff);
  class cov_pg;
  using roffs_t = std::vector<result_offset>;

  friend class read_cov_benchmark;
  static void translate_uint32s(uint32_t* start, uint32_t* limit,
                                const std::vector<uint32_t>& translate_table);

  roff_index_t add_result_offset(result_offset roff);
  void gc_results();
  void advance_to(aoffset_t offset);
  void flush_active_to_here();
  void flush();
  std::unique_ptr<cov_pg> new_pg();
  void on_add_stats();
  void display_stats_report();

  std::unique_ptr<cov_pg> process_assembly(std::unique_ptr<cov_pg> var_pg, assembly_ptr a);

  aoffset_t m_cur_offset = 0;
  bool m_did_notify_path_trim = false;
  assemble_options m_options;
  size_t m_result_idx = 0;

  std::map<aoffset_t /* right offset */, std::unique_ptr<cov_pg>> m_active;
  std::map<aoffset_t /* right offset */, std::unique_ptr<cov_pg>> m_ref_active;

  // Counters of current number of structures allocated
  size_t m_tot_pg = 0;
  size_t m_tot_results = 0;
  size_t m_tot_roffs = 0;
  size_t m_tot_paths = 0;
  size_t m_tot_roff_refs = 0;
  size_t m_tot_pialloc_blocks = 0;
  size_t m_tot_big_pialloc_blocks = 0;
  size_t m_tot_big_pialloc_size = 0;

  size_t m_assemblies_seen = 0;
  time_t m_last_stats_report = 0;

  // Unused blocks of roff indexes, in chunks of k_pi_alloc_chunk elements.
  // 2MB is a TLB entry, so allocate in 2 MB chunks
  static constexpr size_t k_pi_alloc_chunk = 2 * 1024 * 1024 / sizeof(roff_index_t);
  // How big the free chunk list should get per element in m_active
  // before clearing it out.
  static constexpr size_t k_max_pi_free_blocks_per_active = 32;
  std::vector<std::unique_ptr<roff_index_t[]>> m_free_pi_allocs;

  // Current assemblies with reference length 0.  These need to be tracked separately because they
  // will branch off and rejoin at the exact same location.
  std::vector<assembly_ptr> m_cur_inserts;
  // Current assemblies with reference length >0 that will rejoin at a
  // later location.
  std::vector<assembly_ptr> m_cur_non_inserts;

  // Results that have been created that may still be referenced by some assemblies.
  absl::flat_hash_set<std::unique_ptr<result>> m_results;
  time_t m_last_results_gc = 0;
};

}  // namespace variants

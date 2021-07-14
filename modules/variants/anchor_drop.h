#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

class anchor_dropper : public sorted_output_pipeline_step {
 public:
  anchor_dropper(const assemble_options& options, pipeline_step_t output);

  void on_assembly(assembly_ptr a) override;

  static unsigned set_max_kmer_size_for_testing(unsigned kmer_size) {
    auto old = g_max_kmer_size;
    g_max_kmer_size = kmer_size;
    return old;
  }

 private:
  struct kmer_it {
    scaffold::iterator scaffold_it;
    unsigned kmer_needed;
    kmer_t kmer;
  };

  void advance_read_ahead_to(aoffset_t offset);
  void advance_trail_behind_to(aoffset_t offset);
  void advance_range_to(aoffset_t trail_behind, aoffset_t read_ahead);

  void advance_kmer(
      kmer_it& pos, aoffset_t offset,
      const std::function<void(kmer_t, aoffset_t)>& process_kmer_f) const;
  void skip_to(kmer_it& pos, aoffset_t offset) const;
  bool try_long_rejoin(assembly_ptr& a);

  static unsigned g_max_kmer_size;
  assemble_options m_options;
  const scaffold* m_scaffold = nullptr;
  unsigned m_kmer_size;

  kmer_t m_kmer_mask;
  aoffset_t m_kmer_skip_len;

  kmer_it m_read_ahead, m_trail_behind;
  std::unordered_multimap<kmer_t, aoffset_t, unsalted_hash> m_kmers;
};

}  // namespace variants

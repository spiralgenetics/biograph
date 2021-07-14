#pragma once

#include "modules/bio_format/make_vars.h"
#include "modules/variants/assemble.h"

namespace variants {

class call_variants {
 public:
  call_variants(const reference* ref, const seqset* the_seqset,
                const readmap* readmap);

  const manifest& get_manifest() const { return m_var_out; }

  void process_assembly(const assembly& assembly,
                        std::vector<struct_var>* struct_vars = nullptr);

  void flush();

 private:
  struct kvsink_entry {
    std::unique_ptr<kv_sink> sink;
    std::unique_ptr<manifest> local_manifest;
  };

  kvsink_entry acquire_kvsink();
  void release_kvsink(kvsink_entry kvp);

  // Counters
  std::atomic<size_t> m_duplicates{0};
  std::atomic<size_t> m_too_many_anchors{0};
  std::atomic<size_t> m_multiple_all_sv{0};
  std::atomic<size_t> m_passed_filter{0};
  std::atomic<size_t> m_failed_filter{0};

  std::mutex m_mutex;
  std::vector<kvsink_entry> m_kvps;
  std::atomic<size_t> m_manifest_index{0};
  manifest m_var_out;

  const reference* m_ref = nullptr;
  const seqset* m_seqset = nullptr;
  const readmap* m_readmap = nullptr;
  std::set<dna_sequence> m_seen_sequences;
};

}  // namespace variants


#include <unordered_map>
#include "modules/bio_base/kmer.h"
#include "modules/io/track_mem.h"

typedef std::pair<kmer_t, uint32_t> overrep_t;

class overrep_map {
 public:
  overrep_map(size_t kmer_size);
  ~overrep_map();
  void add_overrep(const overrep_t& k);
  bool find_near(const kmer_t& k, overrep_t& out) const;
  size_t size() const { return m_overreps.size(); }
  bool empty() const { return m_overreps.empty(); }

 private:
  void try_side(const kmer_t& k, overrep_t& out) const;
  size_t m_kmer_size;
  tracked_vector<overrep_t> m_overreps;
  tracked_unordered_multimap<uint32_t, uint32_t> m_half0;
  tracked_unordered_multimap<uint32_t, uint32_t> m_half1;
};

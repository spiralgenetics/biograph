#include "modules/bio_base/biograph.h"

biograph::biograph(const std::string& dirname, cache_strategy strategy)
    : m_bgdir(dirname, READ_BGDIR), m_strategy(strategy) {
  if (m_strategy == cache_strategy::RAM) {
    m_options.read_into_ram = true;
  }
}

std::shared_ptr<readmap> biograph::open_readmap(const std::string& id) {
  std::string readmap_path = m_bgdir.find_readmap(id);
  auto& weak = m_readmaps[readmap_path];
  std::shared_ptr<readmap> result = weak.lock();
  if (!result) {
    result = std::make_shared<readmap>(get_seqset(), readmap_path, m_options);
    weak = result;
    if (m_strategy == cache_strategy::MMAPCACHE) {
      result->membufs().cache_in_memory();
    }
    result->calc_read_len_limits_if_needed();
  }
  return result;
}

std::shared_ptr<seqset> biograph::get_seqset() {
  if (!m_seqset) {
    m_seqset.reset(new seqset(m_bgdir.seqset(), m_options));
    if (m_strategy == cache_strategy::MMAPCACHE) {
      m_seqset->membufs().cache_in_memory();
    }
    if (m_strategy != cache_strategy::MMAP) {
      // Often times API users use multiprocessing instead of
      // multithreading, so we want to compute the entry_shared
      // summary table before we fork.  But only do this if we're
      // reading the seqset into RAM anyways.
      m_seqset->init_shared_lt_search();
    }
  }
  return m_seqset;
}

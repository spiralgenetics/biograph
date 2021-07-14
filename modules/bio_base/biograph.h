#pragma once

#include "modules/bio_base/biograph.h"
#include "modules/bio_base/biograph_dir.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"

// A wrapper around biograph_dir that provides seqset and readmap objects.

class biograph {
 public:
  enum class cache_strategy { MMAP, MMAPCACHE, RAM };

  biograph(const std::string& dirname, cache_strategy strategy = cache_strategy::MMAPCACHE);
  ~biograph() = default;

  // Opens a readmap with the given ID, which can either be a UUID, an
  // accession_id, or an empty string.  If an empty string is
  // provided, there must be only one readmap present, and that
  // singular readmap is returned.
  std::shared_ptr<readmap> open_readmap(const std::string& id="");

  // Opens the seqset associated with this biograph.
  std::shared_ptr<seqset> get_seqset();

  const biograph_metadata& get_metadata() const { return m_bgdir.get_metadata(); }

 private:
  biograph_dir m_bgdir;
  cache_strategy m_strategy;
  spiral_file_options m_options;

  std::map<std::string /* path */, std::weak_ptr<readmap>> m_readmaps;
  std::shared_ptr<seqset> m_seqset;
};

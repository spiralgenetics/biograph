#pragma once

#include "modules/io/bitcount.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file.h"
#include "modules/io/version.h"

extern const product_version k_mergemap_version;

struct mergemap_metadata {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(orig_seqset_uuid);
    FIELD(merged_seqset_uuid);
  }

  std::string orig_seqset_uuid;
  std::string merged_seqset_uuid;
};

class seqset_mergemap {
 public:
  seqset_mergemap(const spiral_file_open_state &state);

  const mergemap_metadata metadata() const { return m_metadata; }

  const bitcount &get_bitcount() const { return *m_merged_entries; }

 private:
  mergemap_metadata m_metadata;
  std::unique_ptr<bitcount> m_merged_entries;
};

class seqset_mergemap_builder {
 public:
  seqset_mergemap_builder(const spiral_file_create_state &state,
                          const std::string &orig_seqset_uuid,
                          const std::string &merged_seqset_uuid,
                          size_t entry_count);

  void finalize(progress_handler_t progress = null_progress_handler);

  void set(size_t index) { m_merged_entries->set(index, true); }

 private:
  std::unique_ptr<bitcount> m_merged_entries;
};

#include "modules/bio_base/seqset_mergemap.h"

const product_version k_mergemap_version{"1.0.0"};

seqset_mergemap_builder::seqset_mergemap_builder(
    const spiral_file_create_state& state, const std::string& orig_seqset_uuid,
    const std::string& merged_seqset_uuid, size_t merged_entry_count) {
  state.set_version("mergemap", k_mergemap_version);

  mergemap_metadata metadata;
  metadata.orig_seqset_uuid = orig_seqset_uuid;
  metadata.merged_seqset_uuid = merged_seqset_uuid;
  state.create_json<mergemap_metadata>("mergemap.json", metadata);

  m_merged_entries.reset(
      new bitcount(state.create_subpart("merged_entries"), merged_entry_count));
}

void seqset_mergemap_builder::finalize(progress_handler_t progress) {
  m_merged_entries->finalize(progress);
}

seqset_mergemap::seqset_mergemap(const spiral_file_open_state& state) {
  state.enforce_max_version("mergemap", k_mergemap_version);

  m_metadata = state.open_json<mergemap_metadata>("mergemap.json");

  m_merged_entries.reset(new bitcount(state.open_subpart("merged_entries")));
}

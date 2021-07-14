#pragma once

#include "modules/bio_base/seqset.h"
#include "modules/build_seqset/part_repo.h"
#include "modules/io/track_mem.h"

namespace build_seqset {

class builder {
 public:
  void build_chunks(part_repo& entries, const std::string& pass_name, bool keep_tmp = true,
                    progress_handler_t progress = null_progress_handler);

  std::unique_ptr<seqset> make_seqset(const spiral_file_create_state& state,
                                      progress_handler_t progress = null_progress_handler);

 private:
  struct built_chunk {
    dna_base_array<boost::optional<tracked_vector<bool>>> has_prev;
    boost::optional<mutable_packed_varbit_vector> sizes;
    boost::optional<mutable_packed_varbit_vector> shared;

    built_chunk() {
      for (auto& v : has_prev) {
        v.emplace(track_alloc("built_chunk:has_prev"));
      }
    }
    built_chunk(const built_chunk&) = delete;
  };

  std::map<dna_sequence /* prefix */, builder::built_chunk> m_chunks;
  unsigned m_max_read_len = 0;
};

}  // namespace build_seqset

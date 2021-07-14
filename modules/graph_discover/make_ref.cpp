#include "modules/graph_discover/make_ref.h"
#include "modules/bio_base/kmer.h"

namespace variants {

constexpr char k_make_ref_name[] = "MAKE_REF";

std::vector<assembly_ptr> make_ref_assemblies(const scaffold& s, aoffset_t start_offset,
                                              aoffset_t end_offset, aoffset_t max_chunk_size) {
  std::vector<assembly_ptr> ref_assemblies;
  scaffold sub = s.subscaffold(start_offset, end_offset - start_offset);
  for (const auto& ext : sub.extents()) {
    aoffset_t ext_pos = ext.offset;
    aoffset_t ext_end = ext.offset + ext.sequence.size();

    while (ext_pos < ext_end) {
      aoffset_t chunk_size = ext_end - ext_pos;
      if (max_chunk_size && chunk_size > max_chunk_size) {
        chunk_size = max_chunk_size;
      }

      assembly_ptr a = make_unique<assembly>();
      a->assembly_id = allocate_assembly_id();
      a->left_offset = start_offset + ext_pos;
      a->seq = ext.sequence.subseq(ext_pos - ext.offset, chunk_size);
      a->right_offset = start_offset + ext_pos + a->seq.size();
      a->matches_reference = true;
      a->tags.insert(k_make_ref_name);
      ref_assemblies.push_back(std::move(a));

      ext_pos += chunk_size;
    }
  }
  return ref_assemblies;
}

}  // namespace variants

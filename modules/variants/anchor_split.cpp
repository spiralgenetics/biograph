#include "modules/variants/anchor_split.h"

#include "modules/variants/scaffold.h"

namespace variants {

void anchor_splitter::on_assembly(assembly_ptr a) {
  flush_sorted_to(a->left_offset);
  if (a->matches_reference) {
    sort_and_output(std::move(a));
    return;
  }

  CHECK_LE(a->right_anchor_len + a->left_anchor_len, a->seq.size()) << *a;

  if (a->left_anchor_len) {
    CHECK_LE(a->left_anchor_len, (a->right_offset - a->left_offset));
    CHECK_LE(a->left_anchor_len, a->seq.size());

    aoffset_t split_pos = a->left_anchor_len;
    auto split = split_assembly(std::move(a), split_pos, split_pos);
    if (m_options.trace_reference_assemblies) {
      split.first->matches_reference = true;
      split.first->left_anchor_len = 0;
      split.first->right_anchor_len = 0;
      sort_and_output(std::move(split.first));
    }
    a = std::move(split.second);
    CHECK_EQ(a->left_anchor_len, 0);
  }

  if (a->right_anchor_len) {
    CHECK_LE(a->right_anchor_len, (a->right_offset - a->left_offset));
    CHECK_LE(a->right_anchor_len, a->seq.size());

    aoffset_t seq_split_pos = aoffset_t(a->seq.size()) - a->right_anchor_len;
    aoffset_t ref_split_pos =
        (a->right_offset - a->left_offset) - a->right_anchor_len;
    auto split = split_assembly(std::move(a), seq_split_pos, ref_split_pos);
    if (m_options.trace_reference_assemblies) {
      split.second->matches_reference = true;
      split.second->left_anchor_len = 0;
      split.second->right_anchor_len = 0;
      sort_and_output(std::move(split.second));
    }
    a = std::move(split.first);
    CHECK_EQ(a->right_anchor_len, 0);
  }

  sort_and_output(std::move(a));
}

}  // namespace variants

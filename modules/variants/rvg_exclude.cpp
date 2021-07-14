#include "modules/variants/rvg_exclude.h"

namespace variants {

void rvg_exclude::flush() {
  for (auto& elem : m_backlog) {
    for (auto& a : elem.second) {
      untrack_left_offset(a->left_offset);
    }
  }
  flush_sorted_to(std::numeric_limits<aoffset_t>::max());
}

void rvg_exclude::on_assembly(assembly_ptr a) {
  flush_sorted_to(a->left_offset);

  if (a->matches_reference) {
    sort_and_output(std::move(a));
    return;
  }

  unsigned reflen = a->right_offset - a->left_offset;
  unsigned seqlen = a->seq.size();

  if (!reflen || !seqlen) {
    // Compare as if VCF padded
    ++reflen;
    ++seqlen;
  }

  if (seqlen >= m_options.vcf_sv_size_threshold || reflen >= m_options.vcf_sv_size_threshold) {
    // Structural variant; write this one out.
    m_known_inphase.insert(a->assembly_id);
    auto it = m_backlog.find(a->assembly_id);
    if (it != m_backlog.end()) {
      for (auto& old_a : it->second) {
        untrack_left_offset(old_a->left_offset);
        sort_and_output(std::move(old_a));
      }
      m_backlog.erase(it);
    }
    sort_and_output(std::move(a));
    return;
  }

  bool output_this = false;
  if (a->other_pair_depth || m_known_inphase.count(a->assembly_id)) {
    output_this = true;
  } else {
    CHECK(!a->pair_coverage.empty());
    output_this = true;
    for (int depth : a->pair_coverage) {
      if (!depth) {
        output_this = false;
        break;
      }
    }
  }

  if (output_this) {
    sort_and_output(std::move(a));
    return;
  }

  track_left_offset(a->left_offset);
  m_backlog[a->assembly_id].push_back(std::move(a));
}

}  // namespace variants

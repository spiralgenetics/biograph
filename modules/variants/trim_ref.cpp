#include "modules/variants/trim_ref.h"
#include "modules/bio_base/kmer.h"

#include <boost/range/irange.hpp>

namespace variants {

static constexpr bool k_trim_debug = false;
constexpr aoffset_t ref_trimmer::k_max_backtrack_len;

void ref_trimmer::flush() { flush_sorted(); }

void ref_trimmer::on_assembly(assembly_ptr a) {
  bool both_anchors = a->left_offset && a->right_offset;
  flush_sorted_to(min(a->left_offset, a->right_offset) - k_max_backtrack_len);
  if (a->matches_reference) {
    if (k_trim_debug) {
      std::cout << *a << " matches reference entirely; skipping\n";
    }
    return;
  }

  scaffold s;
  if (both_anchors) {
    s = m_scaffold->subscaffold(a->left_offset, a->right_offset - a->left_offset);
  } else if (a->left_offset) {
    s = m_scaffold->subscaffold(a->left_offset, m_scaffold->end_pos() - a->left_offset);
  } else {
    CHECK(a->right_offset);
    s = m_scaffold->subscaffold(0, a->right_offset);
  }

  if (k_trim_debug) {
    std::cout << "ref_trimmer processing assembly " << *a << " scaffold " << s << "\n";
  }

  if (s.empty()) {
    if (k_trim_debug) {
      std::cout << "ref_trimmer encountered assembly without reference: " << *a << "\n";
    }
    sort_and_output(std::move(a));
    return;
  }
  int shared_left = 0;
  if (a->left_offset) {
    shared_left = s.shared_prefix_length(dna_slice(a->seq));
    CHECK_GE(shared_left, a->left_anchor_len) << *a;
  }

  int shared_right = 0;
  if (a->right_offset) {
    shared_right = s.rev_comp().shared_prefix_length(dna_slice(a->seq).rev_comp());
    CHECK_GE(shared_right, a->right_anchor_len) << *a;
  }

  aoffset_t max_anchor_size = a->seq.size();
  if (both_anchors) {
    max_anchor_size = std::min<aoffset_t>(max_anchor_size, a->right_offset - a->left_offset);
  }

  if (shared_right + shared_left > max_anchor_size) {
    if (shared_right + a->left_anchor_len <= max_anchor_size) {
      shared_left = max_anchor_size - shared_right;
    }
  }
  if (shared_right + shared_left > max_anchor_size) {
    if (a->right_anchor_len + shared_left <= max_anchor_size) {
      shared_right = max_anchor_size - shared_left;
    }
  }

  if (shared_right + shared_left <= max_anchor_size) {
    a->left_anchor_len = shared_left;
    a->right_anchor_len = shared_right;
  }

  if (!a->right_offset && a->right_anchor_len > k_max_backtrack_len) {
    a->right_anchor_len = k_max_backtrack_len;
  }

  aoffset_t new_size = aoffset_t(a->seq.size()) - a->left_anchor_len - a->right_anchor_len;
  CHECK_GE(new_size, 0);

  a->seq = a->seq.subseq(a->left_anchor_len, new_size);
  a->seqset_entries.clear();
  a->rc_seqset_entries.clear();
  if (a->left_offset) {
    a->left_offset += a->left_anchor_len;
    a->left_anchor_len = 0;
  } else {
    CHECK_EQ(a->left_anchor_len, 0);
  }

  if (a->right_offset) {
    a->right_offset -= a->right_anchor_len;
    a->right_anchor_len = 0;
  } else {
    CHECK_EQ(a->right_anchor_len, 0);
  }

  if (both_anchors) {
    CHECK_GE(a->right_offset, a->left_offset);
  }

  if (new_size == 0 && (!both_anchors || a->left_offset == a->right_offset)) {
    if (k_trim_debug) {
      std::cout << "Ref_Trimmer dropping variant that entirely matches reference\n";
    }
    return;
  }

  sort_and_output(std::move(a));
}

}  // namespace variants

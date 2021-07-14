#include "modules/variants/normalize.h"

#include "modules/bio_base/seqset.h"
#include "modules/variants/scaffold.h"

namespace variants {

constexpr bool k_norm_dbg = false;

void normalizer::flush() { flush_sorted_to(std::numeric_limits<aoffset_t>::max()); }

void normalizer::on_assembly(assembly_ptr a) {
  if (a->matches_reference) {
    sort_and_output(std::move(a));
    return;
  }
  a->aligned_variants.clear();
  CHECK(a->aligned_variants.empty());

  if (k_norm_dbg) {
    std::cout << "Normalizer received assembly: " << *a << "\n";
  }

  for (;;) {
    if (k_norm_dbg) {
      std::cout << "Normalizer normalizing: " << *a << "\n";
    }
    auto left_split = m_options.scaffold->split_extent_at(a->left_offset);
    dna_slice left_ref = left_split.first;
    advance_to(a->left_offset - aoffset_t(left_ref.size()) - 1);

    if (k_norm_dbg) {
      dna_slice left_ref_print = left_ref;
      if (left_ref_print.size() > 50) {
        left_ref_print = left_ref_print.subseq(left_ref_print.size() - 50, 50);
      }
      dna_slice right_ref_print = left_split.second;
      if (right_ref_print.size() > 50) {
        right_ref_print = right_ref_print.subseq(0, 50);
      }

      std::cout << "Reference before start: " << left_ref_print << "\n";
      std::cout << "Reference after start: " << left_ref_print << "\n";
    }
    if (left_ref.size() == 0 ||
        aoffset_t(left_split.second.size()) < (a->right_offset - a->left_offset)) {
      if (k_norm_dbg) {
        std::cout << "Too small; outputting " << *a << "\n";
      }
      sort_and_output(std::move(a));
      return;
    }

    dna_slice asm_ref = left_split.second.subseq(0, a->right_offset - a->left_offset);
    dna_base left_base = left_ref[left_ref.size() - 1];
    boost::optional<dna_base> right_asm_base, right_ref_base;
    // Work around gcc bug
    right_asm_base.emplace(int(0));
    right_asm_base.reset();

    if (asm_ref.size() != 0) {
      right_asm_base.emplace(asm_ref[asm_ref.size() - 1]);
    }
    if (a->seq.size() != 0) {
      right_ref_base.emplace(a->seq[a->seq.size() - 1]);
    }
    CHECK(right_asm_base || right_ref_base);
    dna_base right_base;
    if (right_asm_base && right_ref_base) {
      if (*right_asm_base != *right_ref_base) {
        if (k_norm_dbg) {
          std::cout << "Done normalizing; asm base " << *right_asm_base << " != ref base "
                    << *right_ref_base << "; outputting " << *a << "\n";
        }
        sort_and_output(std::move(a));
        return;
      }
      right_base = *right_asm_base;
    }

    if (right_asm_base) {
      right_base = *right_asm_base;
    } else {
      CHECK(right_ref_base);
      right_base = *right_ref_base;
    }

    if (left_base != right_base) {
      if (k_norm_dbg) {
        std::cout << "Done normalizing; left base " << left_base << " !=  right base " << right_base
                  << "; outputting " << *a << "\n";
      }
      sort_and_output(std::move(a));
      return;
    }

    a->left_offset--;
    a->right_offset--;
    if (a->seq.size() > 0) {
      dna_sequence new_seq;
      new_seq.push_back(left_base);
      new_seq += a->seq.subseq(0, a->seq.size() - 1);
      a->seq = new_seq;
      if (k_norm_dbg) {
        std::cout << "Normalized: " << *a << "\n";
      }
    } else {
      if (k_norm_dbg) {
        std::cout << "Moved left: " << *a << "\n";
      }
    }
  }

  // Shift inserts left
  // "TC" + "T" + "A" -> "TC" + "GGC" + "A"
  // turns into:
  // "T" + "CGG" + "CA"
}

void normalizer::advance_to(aoffset_t new_offset) {
  if (new_offset < m_cur_offset) {
    return;
  }
  m_cur_offset = new_offset;
  flush_sorted_to(new_offset);
}

void vcf_padder::on_assembly(assembly_ptr a) {
  if (a->left_offset > 3) {
    flush_sorted_to(a->left_offset - 2);
  }

  if (a->left_offset < a->right_offset && a->seq.size() > 0) {
    sort_and_output(std::move(a));
    return;
  }

  if (a->left_offset <= 0) {
    // Can't do anything here.
    sort_and_output(std::move(a));
    return;
  }

  scaffold s = m_options.scaffold->subscaffold(a->left_offset - 1, 1);
  if (!s.is_simple()) {
    // Can't do anything here.
    sort_and_output(std::move(a));
    return;
  }

  --a->left_offset;
  dna_slice ref_base = s.get_simple();
  CHECK_EQ(ref_base.size(), 1);
  dna_sequence new_seq;
  new_seq.push_back(ref_base[0]);
  new_seq += a->seq;
  a->seq = new_seq;

  sort_and_output(std::move(a));
  return;
}

void vcf_padder::flush() { flush_sorted_to(std::numeric_limits<aoffset_t>::max()); }

}  // namespace variants

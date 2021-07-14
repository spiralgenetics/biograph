#include "modules/variants/dedup.h"

#include "modules/variants/scaffold.h"

namespace variants {

constexpr bool k_dedup_debug = false;

assembly_ptr deduper::combine(assembly_ptr a, assembly_ptr target) {
  if (k_dedup_debug) {
    std::cout << "Combining " << *a << " into: " << *target << "\n";
  }

  // Make sure "a" is to the right of the target; this way we don't
  // have to insert bases or coverage at the beginning of the target.
  if (a->left_offset < target->left_offset) {
    if (k_dedup_debug) {
      std::cout << "Swapping to pad\n";
    }
    std::swap(a, target);
  }

  aoffset_t right_seq_to_add = a->right_offset - target->right_offset;
  if (right_seq_to_add > 0) {
    CHECK_LT(right_seq_to_add, a->seq.size());
    target->seq += a->seq.subseq(a->seq.size() - right_seq_to_add, right_seq_to_add);
    target->right_anchor_len += right_seq_to_add;
    target->right_offset += right_seq_to_add;
  }

  if (!a->coverage.empty() || !target->coverage.empty()) {
    target->coverage.resize(target->seq.size(), 0);

    aoffset_t coverage_offset = a->left_offset - target->left_offset;
    CHECK_GE(coverage_offset, 0);

    if (!a->coverage.empty()) {
      for (unsigned i = 0; i < a->coverage.size(); ++i) {
        target->coverage[i + coverage_offset] += a->coverage[i];
      }
    }
  }

  for (auto matches : {&assembly::left_pair_matches, &assembly::right_pair_matches}) {
    auto& pm = (*target).*matches;
    const auto& old_pm = (*a).*matches;
    pm.insert(pm.end(), old_pm.begin(), old_pm.end());
    std::sort(pm.begin(), pm.end());
    pm.erase(std::unique(pm.begin(), pm.end()), pm.end());
  }

  target->merged_assembly_ids.push_back(a->assembly_id);
  target->merged_assembly_ids.insert(target->merged_assembly_ids.end(),
                                     a->merged_assembly_ids.begin(), a->merged_assembly_ids.end());

  CHECK_GE(target->left_anchor_len, a->left_anchor_len);
  CHECK_GE(target->right_anchor_len, a->right_anchor_len);

#define MERGE_VAR(VARNAME)                                   \
  do {                                                       \
    target->VARNAME = std::max(target->VARNAME, a->VARNAME); \
  } while (0)

  MERGE_VAR(trace_steps);
  MERGE_VAR(unique_pairs_used);
  MERGE_VAR(min_overlap);
  MERGE_VAR(left_anchor_ambiguous_bases);
#undef MERGE_VAR

  if (target->ml_features || a->ml_features) {
    CHECK(target->ml_features);
    CHECK(a->ml_features);

    if (target->ml_features->alt_seq.size() < a->ml_features->alt_seq.size()) {
      target->ml_features = a->ml_features;
    }
  }
  return target;
}

void deduper::on_assembly(assembly_ptr a) {
  if (k_dedup_debug) {
    std::cout << "Deduper input: " << *a << "\n";
  }
  advance_to(a->left_offset);

  auto it = m_queued.begin();
  for (it = m_queued.begin(); it != m_queued.end(); ++it) {
    const assembly_ptr& q = it->second;

    if (k_dedup_debug) {
      std::cout << "Comparing " << *q << " against " << *a << "\n";
    }

    if (q->left_offset > a->left_offset + a->left_anchor_len) {
      // Can't merge any more after this; stop searching.
      break;
    }
    if (a->matches_reference || q->matches_reference) {
      // Don't merge reference-only assemblies.
      continue;
    }

    // Ref trimming should be done already, so only merge if the
    // non-ref portion of the assembly starts and ends in the same
    // place.
    if ((a->left_offset + a->left_anchor_len) != (q->left_offset + q->left_anchor_len)) {
      continue;
    }

    if (a->right_offset - a->right_anchor_len != q->right_offset - q->right_anchor_len) {
      continue;
    }

    CHECK_GE(a->seq.size(), a->left_anchor_len + a->right_anchor_len) << *a;
    dna_slice avar = dna_slice(a->seq).subseq(
        a->left_anchor_len, a->seq.size() - (a->left_anchor_len + a->right_anchor_len));
    CHECK_GE(q->seq.size(), q->left_anchor_len + q->right_anchor_len) << *a;
    dna_slice qvar = dna_slice(q->seq).subseq(
        q->left_anchor_len, q->seq.size() - (q->left_anchor_len + q->right_anchor_len));
    if (avar != qvar) {
      continue;
    }
    if (k_dedup_debug) {
      std::cout << "Compared assemblies: " << *a << " vs " << *q << "\n";
      std::cout << "Compared variants: " << avar << " vs " << qvar << "\n";
      std::cout << "Can merge!\n";
    }

    a = combine(std::move(it->second), std::move(a));
    m_queued.erase(it);
    break;
  }

  auto left_offset = a->left_offset;
  m_queued.emplace(left_offset, std::move(a));
}

void deduper::advance_to(aoffset_t offset) {
  while (!m_queued.empty()) {
    auto it = m_queued.begin();
    if (it->second->left_offset + it->second->left_anchor_len >= offset) {
      // Still a candidate for merging; leave the rest in queue.
      break;
    }
    m_output->add(std::move(it->second));
    m_queued.erase(it);
  }
}

void deduper::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_queued.empty());
}

void exact_deduper::on_assembly(assembly_ptr a) {
  if (k_dedup_debug) {
    std::cout << "Exact deduper input: " << dump_assembly_and_vars(*a) << "\n";
  }
  advance_to(a->left_offset);

  auto it = m_queued.begin();
  for (it = m_queued.begin(); it != m_queued.end(); ++it) {
    assembly_ptr& q = it->second;

    if (a->left_offset == q->left_offset && a->right_offset == q->right_offset &&
        a->seq == q->seq && a->matches_reference == q->matches_reference) {
      if (k_dedup_debug) {
        std::cout << "Mostly discarding: " << *a << "\nIn favor of: " << *q << "\n";
      }
      for (auto matches : {&assembly::left_pair_matches, &assembly::right_pair_matches}) {
        auto& pm = (*q).*matches;
        const auto& old_pm = (*a).*matches;
        pm.insert(pm.end(), old_pm.begin(), old_pm.end());
        std::sort(pm.begin(), pm.end());
        pm.erase(std::unique(pm.begin(), pm.end()), pm.end());
      }

      q->merged_assembly_ids.push_back(a->assembly_id);
      for (const auto& id : a->merged_assembly_ids) {
        q->merged_assembly_ids.push_back(id);
      }

      if (q->ml_features || a->ml_features) {
        CHECK(q->ml_features);
        CHECK(a->ml_features);
        if (q->ml_features->alt_seq.size() < a->ml_features->alt_seq.size()) {
          q->ml_features = a->ml_features;
        }
      }

      if (k_dedup_debug) {
        std::cout << "After merge: " << *q << "\n";
      }
      return;
    }
  }

  auto left_offset = a->left_offset;
  m_queued.emplace(left_offset, std::move(a));
}

void exact_deduper::advance_to(aoffset_t offset) {
  while (!m_queued.empty()) {
    auto it = m_queued.begin();
    if (it->first >= offset) {
      return;
    }

    if (k_dedup_debug) {
      std::cout << "Exact dedupper outputting " << *it->second << " at " << offset << "\n";
    }
    m_output->add(std::move(it->second));
    m_queued.erase(it);
  }
}

void exact_deduper::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_queued.empty());
}

}  // namespace variants

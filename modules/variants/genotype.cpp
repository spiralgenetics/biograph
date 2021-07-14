#include "modules/variants/genotype.h"

namespace variants {

constexpr bool k_gt_dbg = false;

std::ostream& operator<<(std::ostream& os, const genotyper::entry& e) {
  if (e.active) {
    os << "ACTIVE: ";
  } else {
    os << "inactive: ";
  }
  os << "DP:" << e.cur_depth << " seq=" << e.seq_offset << " ref=" << e.ref_offset;
  if (e.a) {
    os << " " << dump_assembly_and_vars(*e.a);
  } else {
    os << " (null assembly)";
  }
  return os;
}

genotyper::genotyper(const assemble_options& options, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_options(options) {}

void genotyper::flush() {
  m_intake_offset = std::numeric_limits<aoffset_t>::max();
  advance();
  CHECK(m_entries.empty());
}

void genotyper::advance() {
  if (k_gt_dbg) {
    std::cout << "Advancing from " << m_process_offset << " to " << m_intake_offset << "\n";
  }
  while (m_process_offset + 1 < m_intake_offset) {
    advance_base();
  }
  if (k_gt_dbg) {
    std::cout << "Flushing sorted to " << m_intake_offset << "\n";
  }
  if (m_process_offset > 0) {
    flush_sorted_to(m_process_offset);
  }
}

void genotyper::on_assembly(assembly_ptr a) {
  if (k_gt_dbg) {
    std::cout << "Genotyper received assembly " << dump_assembly_and_vars(*a) << "\n";
  }
  CHECK(!a->coverage.empty()) << "Genotyper requires coverage\n";
  CHECK_GE(a->left_offset, m_intake_offset);
  track_left_offset(a->left_offset);
  m_intake_offset = a->left_offset;
  advance();

  entry_ptr e = make_unique<entry>();
  e->a = std::move(a);
  init_entry(*e);
  m_entries.push_back(std::move(e));
}

void genotyper::advance_base() {
  if (m_entries.empty()) {
    // Nothing going on; skip ahead.
    if (k_gt_dbg) {
      std::cout << "Skipping ahead from " << m_process_offset << " to " << m_intake_offset << "\n";
    }
    m_process_offset = m_intake_offset - 1;
    return;
  }

  CHECK_LT(m_process_offset, m_intake_offset);

  if (m_entries.empty()) {
    // Nothing going on; skip ahead.
    if (k_gt_dbg) {
      std::cout << "Skipping ahead from " << m_process_offset << " to " << m_intake_offset << "\n";
    }
    m_process_offset = m_intake_offset - 1;
    return;
  }

  if (k_gt_dbg) {
    std::cout << "Processing at " << m_process_offset << " with " << m_entries.size()
              << " active\n";
  }

  std::stable_sort(
      m_entries.begin(), m_entries.end(),
      [&](const entry_ptr& lhs, const entry_ptr& rhs) { return lhs->cur_depth > rhs->cur_depth; });

  unsigned num_output = 0;
  int best_depth = m_entries.front()->cur_depth;
  for (size_t idx = 0; idx < m_entries.size(); ++idx) {
    const auto& e = m_entries[idx];
    CHECK(e->a);
    CHECK_LE(e->cur_depth, best_depth);
    if (e->cur_depth == 0) {
      deactivate(*e);
      continue;
    }
    if (e->cur_depth < best_depth * m_options.min_depth_portion) {
      deactivate(*e);
      continue;
    }

    if (num_output >= m_options.max_ploids) {
      deactivate(*e);
      continue;
    }

    bool is_ref = e->vit == e->a->aligned_variants.end();
    bool is_duplicate = false;
    for (size_t idx2 = 0; idx2 < idx; ++idx2) {
      const auto& e2 = m_entries[idx2];

      bool is_ref2 = e2->vit == e2->a->aligned_variants.end();

      if (is_ref && is_ref2) {
        is_duplicate = true;
      } else if (!is_ref && !is_ref2) {
        is_duplicate = *e->vit == *e2->vit;
      }

      if (is_duplicate) {
        break;
      }
    }

    if (is_duplicate) {
      deactivate(*e);
      continue;
    }

    if (!e->active && e->in_variant) {
      // Don't activate in the middle of a variant.
      continue;
    }

    if (k_gt_dbg) {
      std::cout << "Active here: id=" << e->a->assembly_id << " depth=" << e->cur_depth << "\n";
    }
    ++num_output;
    activate(*e);
  }

  // Calculate alternate depths
  for (size_t idx = 0; idx < m_entries.size(); ++idx) {
    const auto& e = m_entries[idx];
    if (!e->active || !e->a) {
      continue;
    }
    if (e->vit == e->a->aligned_variants.end() || e->vit->left_offset > m_process_offset) {
      continue;
    }
    int max_alt_depth = 0;
    for (size_t idx2 = 0; idx2 < m_entries.size(); ++idx2) {
      const auto& e2 = m_entries[idx2];
      if (!e2->active || !e2->a || idx == idx2) {
        continue;
      }
      max_alt_depth += m_entries[idx2]->cur_depth;
    }
    if (max_alt_depth > e->vit->max_alt_depth) {
      e->vit->max_alt_depth = max_alt_depth;
    }
  }

  ++m_process_offset;
  size_t idx = 0;
  while (idx < m_entries.size()) {
    advance_entry(*m_entries[idx]);
    if (m_entries[idx]->a) {
      calc_depth(*m_entries[idx]);
      ++idx;
    } else {
      m_entries.erase(m_entries.begin() + idx);
    }
  }
}

bool genotyper::check_entry_done(entry& e) {
  if (e.ref_offset == e.a->right_offset && e.vit == e.a->aligned_variants.end()) {
    if (k_gt_dbg) {
      std::cout << "Entry detected as done: " << *e.a << "\n";
    }
    CHECK_EQ(e.seq_offset, e.a->seq.size());
    assembly_ptr a = std::move(e.a);
    untrack_left_offset(a->left_offset);
    if (is_degenerate(*a)) {
      a.release_and_discard();
    } else {
      if (e.active) {
        sort_and_output(std::move(a));
      } else {
        report_discard(std::move(a));
      }
    }
    return true;
  }
  return false;
}

void genotyper::advance_entry(entry& e) {
  if (!e.a) {
    return;
  }
  if (k_gt_dbg) {
    std::cout << "Advancing entry " << dump_assembly_and_vars(*e.a) << " from seq=" << e.seq_offset
              << ", ref=" << e.ref_offset << " to " << m_process_offset << "\n";
    if (e.vit != e.a->aligned_variants.end()) {
      std::cout << "vit = " << e.vit->left_offset << " to " << e.vit->right_offset << "\n";
    }
  }

  e.in_variant = false;

  if (e.ref_offset < e.a->left_offset) {
    if (k_gt_dbg) {
      std::cout << "Before assembly start\n";
    }
    ++e.ref_offset;
    CHECK_EQ(e.ref_offset, m_process_offset);
    return;
  }

  // If we were in a variant, advance the sequence past it.
  if (e.vit != e.a->aligned_variants.end() && e.vit->right_offset == e.ref_offset) {
    e.seq_offset += e.vit->seq.size();
    ++e.vit;
  }

  if (k_gt_dbg) {
    std::cout << "Before advance, ref=" << e.ref_offset << " seq=" << e.seq_offset << "\n";
  }

  if (check_entry_done(e)) {
    return;
  }

  if (e.vit == e.a->aligned_variants.end() || e.ref_offset < e.vit->left_offset) {
    e.deactivate_seq_offset = e.seq_offset;
    e.deactivate_ref_offset = e.ref_offset;
    ++e.seq_offset;
  } else {
    if (e.ref_offset < e.vit->right_offset) {
      e.in_variant = true;
    }
  }
  ++e.ref_offset;
  CHECK_EQ(e.ref_offset, m_process_offset);

  if (k_gt_dbg) {
    std::cout << "After advance, ref=" << e.ref_offset << " seq=" << e.seq_offset << "\n";
  }

  CHECK_LE(m_process_offset, e.a->right_offset + 1);
}

void genotyper::calc_depth(entry& e) {
  if (m_process_offset < e.a->left_offset || m_process_offset > e.a->right_offset) {
    if (k_gt_dbg) {
      std::cout << "Before first or after last; 0 depth\n";
    }
    e.cur_depth = 0;
    if (m_process_offset < e.a->left_offset) {
      CHECK_EQ(0, e.seq_offset);
    }
    return;
  }

  CHECK_EQ(e.ref_offset, m_process_offset);

  if (e.vit != e.a->aligned_variants.end() && e.vit->left_offset == e.ref_offset) {
    // Calculate depth for the whole variant.
    CHECK(e.vit != e.a->aligned_variants.end());
    int min_depth = std::numeric_limits<int>::max();
    for (aoffset_t idx = 0; idx <= aoffset_t(e.vit->seq.size()); ++idx) {
      int depth = e.a->coverage.at(e.seq_offset + idx);
      if (depth < min_depth) {
        min_depth = depth;
      }
    }
    if (k_gt_dbg) {
      std::cout << "Min depth for variant " << *e.vit << " calculated to be " << min_depth
                << " from " << *e.a << " seq_offset = " << e.seq_offset << "\n";
    }
    e.variant_depth = min_depth;
  }

  if (e.vit != e.a->aligned_variants.end() && m_process_offset >= e.vit->left_offset) {
    CHECK_LE(m_process_offset, e.vit->right_offset);
    e.cur_depth = e.variant_depth;
    if (k_gt_dbg) {
      std::cout << "Depth from variant coverage: " << e.cur_depth << "\n";
    }
  } else {
    CHECK_EQ(m_process_offset, e.ref_offset);
    e.cur_depth = e.a->coverage.at(e.seq_offset);
    if (k_gt_dbg) {
      std::cout << "Depth from id=" << e.a->assembly_id << " non-variant coverage at "
                << e.seq_offset << ": " << e.cur_depth << "\n";
    }
  }
}

void genotyper::report_discard(assembly_ptr a) {
  if (m_options.report_genotype_discard_func) {
    std::vector<const assembly*> active;

    for (const auto& e : m_entries) {
      if (!e->a) {
        continue;
      }
      if (!e->active) {
        continue;
      }
      active.push_back(e->a.get());
    }
    m_options.report_genotype_discard_func(m_options, *a, active);
  }
}

void genotyper::activate(entry& e) {
  if (e.active) {
    return;
  }

  CHECK_EQ(m_process_offset, e.ref_offset);
  aoffset_t ref_split_pos = m_process_offset;
  if (k_gt_dbg) {
    std::cout << "Activating " << dump_assembly_and_vars(*e.a) << " at " << ref_split_pos << "\n";
  }
  aoffset_t seq_split_pos = e.seq_offset;

  if (e.vit != e.a->aligned_variants.end() && ref_split_pos < e.vit->right_offset &&
      ref_split_pos > e.vit->left_offset) {
    // Don't activate in the middle of a variant.
    ref_split_pos = e.vit->right_offset;
    seq_split_pos += e.vit->seq.size();
    ++e.vit;
  }

  aoffset_t left_offset = e.a->left_offset;
  untrack_left_offset(left_offset);
  if (k_gt_dbg) {
    std::cout << "Splitting at seq pos " << seq_split_pos << " ref pos = " << ref_split_pos << "\n";
  }
  aoffset_t ref_len = ref_split_pos - e.a->left_offset;
  auto split = split_assembly(std::move(e.a), seq_split_pos, ref_len);
  if (k_gt_dbg) {
    std::cout << "Split results were " << dump_assembly_and_vars(*split.first) << " and "
              << dump_assembly_and_vars(*split.second) << "\n";
  }
  report_discard(std::move(split.first));
  // Discard split.first; it was not useful.
  e.a = std::move(split.second);
  track_left_offset(e.a->left_offset);
  init_entry(e);
  e.active = true;
}

void genotyper::deactivate(entry& e) {
  if (!e.active) {
    return;
  }

  aoffset_t ref_split_pos = e.deactivate_ref_offset;
  aoffset_t seq_split_pos = e.deactivate_seq_offset;

  if (k_gt_dbg) {
    std::cout << "Deactivating " << dump_assembly_and_vars(*e.a) << " at ref=" << ref_split_pos
              << " seq=" << seq_split_pos << "\n";
  }

  aoffset_t left_offset = e.a->left_offset;
  untrack_left_offset(left_offset);
  if (k_gt_dbg) {
    std::cout << "Splitting at seq pos " << seq_split_pos << " ref pos = " << ref_split_pos << "\n";
  }
  auto ref_len = ref_split_pos - e.a->left_offset;
  auto split = split_assembly(std::move(e.a), seq_split_pos, ref_len);
  if (k_gt_dbg) {
    std::cout << "Split results were " << dump_assembly_and_vars(*split.first) << " and "
              << dump_assembly_and_vars(*split.second) << "\n";
  }
  if (!is_degenerate(*split.first)) {
    sort_and_output(std::move(split.first));
  }
  e.a = std::move(split.second);
  track_left_offset(e.a->left_offset);

  // Roll back and catch this one up past the part we don't want to include.
  aoffset_t orig_process_offset = m_process_offset;
  if (e.a->left_offset < orig_process_offset) {
    m_process_offset = e.a->left_offset;
  }
  init_entry(e);
  while (m_process_offset < orig_process_offset) {
    ++m_process_offset;
    advance_entry(e);
  }
}

void genotyper::init_entry(entry& e) {
  e.active = false;
  e.vit = e.a->aligned_variants.begin();
  e.deactivate_seq_offset = e.seq_offset = 0;
  e.cur_depth = 0;
  e.in_variant = false;
  e.deactivate_ref_offset = e.ref_offset = m_process_offset;
  CHECK_GE(e.a->left_offset, m_process_offset) << *e.a;
  calc_depth(e);
}

bool genotyper::is_degenerate(const assembly& a) {
  return a.seq.size() == 0 && a.left_offset == a.right_offset;
}

}  // namespace variants

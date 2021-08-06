#include "modules/variants/apply_graph.h"

static constexpr bool k_dbg = false;

namespace variants {

apply_graph::apply_graph(const on_context_func_t& f, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_on_context(f) {}

void apply_graph::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
  CHECK(!m_left_ref);
  CHECK(!m_right_ref);
}

void apply_graph::on_assembly(assembly_ptr a) {
  aoffset_t left_offset = min(a->left_offset, a->right_offset);
  aoffset_t right_offset = max(a->left_offset, a->right_offset);
  advance_to(left_offset);

  track_left_offset(left_offset);

  if (a->matches_reference) {
    m_right_ref = aptr::make_shared(std::move(a));
  } else {
    result r;
    r.a = aptr::make_shared(std::move(a));
    if (m_left_ref) {
      r.left_ref = m_left_ref.clone();
    }
    m_active.emplace(right_offset, std::move(r));
  }
}

void apply_graph::advance_to(aoffset_t target) {
  CHECK_GE(target, m_cur_offset);
  while (m_cur_offset < target) {
    advance_towards(target);
  }
}

void apply_graph::advance_towards(aoffset_t target) {
  CHECK_GT(target, m_cur_offset);

  if (k_dbg) {
    std::cerr << "Flushing at " << m_cur_offset << ", advancing towards " << target << "\n";
  }

  // Flush assemblies we're done processing.
  auto act_it = m_active.begin();
  while (act_it != m_active.end()) {
    aoffset_t act_right = act_it->first;
    auto& r = act_it->second;

    CHECK_GE(act_right, m_cur_offset);

    if (act_right > m_cur_offset) {
      break;
    }

    CHECK_EQ(act_right, m_cur_offset);
    // Done with this variant.
    if (m_right_ref) {
      r.right_ref = m_right_ref.clone();
    }
    output_result(std::move(r));
    act_it = m_active.erase(act_it);
  }

  if (!m_active.empty() && m_active.begin()->first < target) {
    target = m_active.begin()->first;
    if (k_dbg) {
      std::cerr << "Target at first active, " << target << "\n";
    }
    CHECK_GT(target, m_cur_offset);
  }

  if (m_right_ref && (*m_right_ref)->right_offset < target) {
    target = (*m_right_ref)->right_offset;
    if (k_dbg) {
      std::cerr << "Target at end of ref, " << target << "\n";
    }
    CHECK_GT(target, m_cur_offset);
  }

  // Fill in reference coverage for any active assemblies.
  for (auto& act_entry : m_active) {
    aoffset_t act_right = act_entry.first;
    auto& r = act_entry.second;

    CHECK_GT(act_right, m_cur_offset);

    aoffset_t act_right_so_far = (*r.a)->left_offset;
    if (!r.refs.empty()) {
      act_right_so_far = (*r.refs.back())->right_offset;
    }

    CHECK_LE(act_right_so_far, m_cur_offset);

    if (m_right_ref) {
      CHECK_LE((*m_right_ref)->right_offset, act_right);

      // Add this reference assembly to the saved refs.
      r.refs.emplace_back(m_right_ref.clone());
    }
  }

  if (m_left_ref) {
    release(std::move(m_left_ref));
  }
  if (m_right_ref) {
    CHECK_EQ((*m_right_ref)->left_offset, m_cur_offset);
    if (k_dbg) {
      std::cerr << "cur=" << m_cur_offset << " target=" << target << " moving ref over? "
                << **m_right_ref << "\n";
    }
    if (target == (*m_right_ref)->right_offset) {
      if (k_dbg) {
        std::cerr << "Yep\n";
      }
      m_left_ref = std::move(m_right_ref);
    } else {
      if (k_dbg) {
        std::cerr << "Nope\n";
      }
      release(std::move(m_right_ref));
    }
  }
  m_cur_offset = target;
}

read_coverage_t graph_context::merge_coverage(
    boost::optional<read_coverage_t> assembly::*field) const {
  CHECK(a);

  if (k_dbg) {
    std::cerr << "Merging for " << *a << "\n";
  }

  CHECK_LE(a->left_offset, a->right_offset);
  aoffset_t reflen = a->right_offset - a->left_offset;

  // Reads that we've verified are contiguous.
  std::vector<read_coverage_read_t> done_reads;

  // Reads that still need to be verified to make sure they continue to the next section
  absl::btree_map<std::pair<aoffset_t /* right offset of read relative to a->left_offset */,
                            aoffset_t /* read len */>,
                  read_coverage_read_t>
      pending_reads, new_pending_reads;

  aoffset_t cur_offset = std::numeric_limits<aoffset_t>::min();

  std::vector<assembly*> outer_refs;
  if (left_ref) {
    if (k_dbg) {
      std::cerr << "Left ref: " << *left_ref << "\n";
    }
    outer_refs.push_back(left_ref);
  }
  outer_refs.insert(outer_refs.end(), refs.begin(), refs.end());
  if (right_ref) {
    if (k_dbg) {
      std::cerr << "Right ref: " << *right_ref << "\n";
    }
    outer_refs.push_back(right_ref);
  }

  for (assembly* ref : outer_refs) {
    CHECK(ref->matches_reference);

    aoffset_t ref_offset = ref->left_offset - a->left_offset;
    CHECK_GE(ref_offset, cur_offset);

    if (k_dbg) {
      std::cerr << "Adding ref " << *ref << " at offset " << ref_offset << " cur=" << cur_offset
                << " with " << pending_reads.size() << " pending and " << done_reads.size()
                << " done\n";
    }

    if (ref_offset > cur_offset) {
      if (k_dbg) {
        std::cerr << "Clearing pending, since ref  is past cur\n";
      }
      pending_reads.clear();
      cur_offset = ref_offset;
    }

    for (const auto& cov_entry : (ref->*field)->reads()) {
      aoffset_t cov_left_offset = cov_entry.offset + cur_offset;
      aoffset_t cov_right_offset = cov_left_offset + cov_entry.read_len;
      if (k_dbg) {
        std::cerr << "Considering entry " << cov_left_offset << " to " << cov_right_offset
                  << ", reflen= " << reflen << ", cur offset = " << cur_offset << "\n";
      }
      if (cov_left_offset >= reflen) {
        continue;
      }
      if (cov_right_offset <= 0) {
        continue;
      }
      std::pair<aoffset_t, aoffset_t> cov_pair{cov_left_offset, cov_entry.read_len};

      read_id_set new_ids;
      if (cov_left_offset < cur_offset && ref_offset >= 0) {
        if (k_dbg){
          std::cerr << "Checking present in " << pending_reads.size() << " prevs\n";}
        // Make sure it's present in the previous reference section.
        auto pending_it = pending_reads.find(cov_pair);
        if (pending_it == pending_reads.end()) {
          // Read not present in previous reference
          if (k_dbg) {
            std::cerr << "not found\n";
          }
          continue;
        }

        new_ids = pending_it->second.read_ids.intersection(cov_entry.read_ids);
        if (new_ids.empty()) {
          if (k_dbg) {
            std::cerr << "no intersection\n";
          }
          continue;
        }
      } else {
        new_ids = cov_entry.read_ids;
      }
      if (k_dbg) {
        std::cerr << new_ids.size() << " reads found\n";
      }

      read_coverage_read_t new_cov_entry;
      new_cov_entry.offset = cov_left_offset;
      new_cov_entry.read_len = cov_entry.read_len;
      new_cov_entry.read_ids = std::move(new_ids);
      if (cov_right_offset <= (ref->right_offset - a->left_offset) ||
          ref->right_offset > a->right_offset) {
        // Read ends in this ref section
        done_reads.emplace_back(std::move(new_cov_entry));
      } else {
        CHECK(new_pending_reads.emplace(cov_pair, std::move(new_cov_entry)).second)
            << "Duplicate coverage position?";
      }
    }

    pending_reads.clear();
    std::swap(pending_reads, new_pending_reads);
    cur_offset = ref->right_offset - a->left_offset;
  }

  if (k_dbg) {
    std::cerr << "Done, " << done_reads.size() << " done and " << pending_reads.size()
              << " still pending\n";
  }

  return read_coverage_t(reflen, std::move(done_reads));
}

void apply_graph::output_result(result r) {
  // Expose context to user
  graph_context ctx;

  ctx.a = r.a->get();

  if (r.left_ref) {
    ctx.left_ref = r.left_ref->get();
  }
  if (r.right_ref) {
    ctx.right_ref = r.right_ref->get();
  }
  for (const auto& ref : r.refs) {
    ctx.refs.push_back(ref->get());
  }

  if (k_dbg) {
    std::cerr << "Outputting result for " << *ctx.a << "\n";
    if (ctx.left_ref) {
      std::cerr << " Left:  " << *ctx.left_ref << "\n";
    }
    for (const auto* ref : ctx.refs) {
      std::cerr << " Mid:   " << ref << "\n";
    }
    if (ctx.right_ref) {
      std::cerr << " Right: " << *ctx.left_ref << "\n";
    }
  }

  // Show the user the context we constructed.
  try {
    m_on_context(ctx);
  } catch (...) {
    clear_result(std::move(r));
    throw;
  }
  clear_result(std::move(r));
}

void apply_graph::clear_result(result r) {
  release(std::move(r.a));
  for (auto& ref : r.refs) {
    release(std::move(ref));
  }
  if (r.left_ref) {
    release(std::move(r.left_ref));
  }
  if (r.right_ref) {
    release(std::move(r.right_ref));
  }
}

void apply_graph::release(aptr a) {
  std::unique_ptr<assembly_ptr> out = a.release();

  if (out) {
    // Last reference; output it.
    untrack_left_offset((*out)->left_offset);
    sort_and_output(std::move(*out));
  }
}

scaffold graph_context::ref_scaffold() const {
  std::vector<scaffold::extent> exts;

  for (const auto* ref : refs) {
    scaffold::extent new_ext;
    new_ext.offset = ref->left_offset - a->left_offset;
    new_ext.sequence = ref->seq;
    exts.push_back(new_ext);
  }

  return scaffold(exts, a->right_offset - a->left_offset);
}

edge_coverage_t graph_context::edge_coverage(const scaffold& s, const read_coverage_t& var_cov,
                                             const read_coverage_t& ref_cov) const {
  edge_coverage_t result;
  aoffset_t reflen = s.end_pos();
  aoffset_t seqlen = a->seq.size();

  // Find how many bases are common with reference at each end.
  result.start_common = s.shared_prefix_length(a->seq);
  result.end_common = s.rev_comp().shared_prefix_length(dna_slice(a->seq).rev_comp());

  // Get variant read counts.
  result.variant_start = var_cov.get_reads_spanning_offset(result.start_common).all_read_ids();
  result.variant_end = var_cov.get_reads_spanning_offset(seqlen - result.end_common).all_read_ids();

  result.reference_start = ref_cov.get_reads_spanning_offset(result.start_common).all_read_ids();
  result.reference_end =
      ref_cov.get_reads_spanning_offset(reflen - result.end_common).all_read_ids();

  // Tally up everything else in the variant into the interior.
  for (const auto& cov_entry : var_cov.reads()) {
    if (cov_entry.offset >= 0 && cov_entry.offset + cov_entry.read_len <= seqlen) {
      result.interior |= cov_entry.read_ids;
    }
  }

  return result;
}

}  // namespace variants

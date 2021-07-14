#include "modules/variants/calc_coverage.h"

#include "modules/bio_base/readmap.h"

namespace variants {

constexpr bool k_cov_dbg = false;

int calc_coverage::sum_dir_depth(int fwd, int rev) const {
  if (m_options.penalize_directional_coverage) {
    int result = fwd + rev;
    int disparity = abs(fwd - rev);
    // If there's a huge disparity between the two directions, discount it.
    if (disparity * 4 > result * 3) {
      return std::min(fwd, rev) * 2;
    }
    return result;
  } else {
    return fwd + rev;
  }
}

bool calc_coverage::coverage_accum::add(int start, int end) {
  CHECK(m_starts);
  CHECK_GT(size(), 0);
  if (start < 0) {
    start = 0;
  }
  if (end >= int(size())) {
    end = size() - 1;
  }
  if (start > int(size() - 1) || end < 0) {
    if (k_cov_dbg) {
      std::cout << "Adjusted start=" << start << " end=" << end << " is out of bounds; skipping\n";
    }
    return false;
  }
  if (k_cov_dbg) {
    std::cout << "Adjusted coverage add: start=" << start << " end=" << end << "\n";
  }
  ++m_starts[start];
  ++m_ends[end];
  return true;
}

std::vector<int> calc_coverage::coverage_accum::coverage() {
  CHECK(m_starts);
  CHECK_GT(size(), 0);
  std::vector<int> result;
  result.reserve(size());
  int cur = 0;
  for (unsigned idx = 0; idx != size() + 1; ++idx) {
    if (k_cov_dbg) {
      std::cout << " ";
    }
    bool changed = false;
    if (idx < size() && m_starts[idx]) {
      if (k_cov_dbg) {
        std::cout << "+" << m_starts[idx];
      }
      cur += m_starts[idx];
      changed = true;
    }
    if (idx > 0 && m_ends[idx - 1]) {
      if (k_cov_dbg) {
        std::cout << "-" << m_ends[idx - 1];
      }
      cur -= m_ends[idx - 1];
      changed = true;
    }
    if (k_cov_dbg) {
      if (changed) {
        std::cout << "=" << cur;
      } else {
        std::cout << cur;
      }
    }
    if (idx == size()) {
      CHECK_EQ(0, cur);
    } else {
      result.push_back(cur);
    }
  }
  if (k_cov_dbg) {
    std::cout << "\n";
  }
  return result;
}

calc_coverage::calc_coverage(const assemble_options& options, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output), true /* old sort order */),
      m_options(options) {
  CHECK(m_options.scaffold);
  CHECK(m_options.readmap);
  set_expected_order(assembly::left_offset_less_than);
  m_scaffold_it = m_options.scaffold->begin();
  m_cur_offset = m_scaffold_it.offset();
  init_ref_pg();
}

calc_coverage::shared_assembly::~shared_assembly() {
  //  std::cout << "Shared assembly destroying: " << *a << "\n";
  if (k_cov_dbg) {
    std::cout << "Shared assembly finishing: " << dump_assembly_and_vars(*a)
              << " at cur_offset=" << c->m_cur_offset << "\n";
    std::cout << "Other mins: " << other_min.size() << "\n";
    for (const auto& cov : other_ref_coverage) {
      std::cout << "Other ref coverage: " << dump_coverage(*cov) << "\n";
    }
  }
  if (other_min.empty()) {
    a->other_depth = 0;
    a->other_pair_depth = 0;
  } else {
    auto it = other_min.begin();
    a->other_depth = it->base;
    a->other_pair_depth = it->pair;
    ++it;

    while (it != other_min.end()) {
      a->other_depth = std::min<int>(a->other_depth, it->base);
      a->other_pair_depth = std::min<int>(a->other_pair_depth, it->pair);
      ++it;
    }
  }
  aoffset_t ref_bases = 0;
  int min_ref_depth = std::numeric_limits<int>::max();
  for (const auto& ref_cov : other_ref_coverage) {
    CHECK(!ref_cov->empty());
    ref_bases += ref_cov->size() - 1;
    min_ref_depth = std::min<int>(min_ref_depth, container_min(*ref_cov));
    if (k_cov_dbg) {
      std::cout << "Got a ref section with coverage length " << ref_cov->size() << "\n";
    }
  }
#ifdef NDEBUG
  constexpr bool k_check_ref_bases = k_cov_dbg ? true : false;
#else
  constexpr bool k_check_ref_bases = true;
#endif
  if (k_check_ref_bases && ref_bases != (a->right_offset - a->left_offset)) {
    if (k_cov_dbg) {
      std::cout << "REF BASE COUNT MISMATCH; got " << ref_bases
                << " bases of ref coverage for assembly " << *a << "\n";
    }
    auto it = c->m_options.scaffold->begin();
    if (it != c->m_options.scaffold->end()) {
      it.skip_to(a->left_offset - 1, "calc_coverage(2)");
      if (it != c->m_options.scaffold->end()) {
        aoffset_t left_extent_end = it.extent_end_offset();
        it.skip_to(a->right_offset, "calc_coverage(3)");
        if (it != c->m_options.scaffold->end()) {
          aoffset_t right_extent_end = it.extent_end_offset();
          CHECK_NE(left_extent_end, right_extent_end);
        }
      }
    }
  }
  CHECK_NE(min_ref_depth, std::numeric_limits<int>::max());
  a->ref_depth = min_ref_depth;
  c->untrack_left_offset(a->left_offset);
  c->sort_and_output(std::move(a));
}

calc_coverage::cov_tracker::~cov_tracker() {
  std::shared_ptr<std::vector<int>> both_dir_shared;
  std::vector<int> both_dir_unshared;
  std::vector<int>* both_dir_cov;
  if (var_assembly) {
    both_dir_cov = &both_dir_unshared;
  } else {
    both_dir_shared = std::make_shared<std::vector<int>>();
    both_dir_cov = both_dir_shared.get();
  }
  auto fwd_cov = fwd_coverage.coverage();
  auto rev_cov = rev_coverage.coverage();
  both_dir_cov->reserve(fwd_cov.size());
  for (unsigned i = 0; i != fwd_cov.size(); ++i) {
    both_dir_cov->push_back(c->sum_dir_depth(fwd_cov[i], rev_cov[i]));
  }

  if (var_assembly) {
    auto pair_cov = pair_coverage.coverage();

    if (!other_assemblies.empty()) {
      other_min_coverage other;
      other.base = container_min(*both_dir_cov);
      other.pair = container_min(pair_cov);

      for (const auto& other_a : other_assemblies) {
        other_a->other_min.push_back(other);
      }
    }

    var_assembly->a->coverage = std::move(*both_dir_cov);
    var_assembly->a->pair_coverage = std::move(pair_cov);
  } else {
    CHECK(!other_assemblies.empty()) << "Useless ref coverage?";
    CHECK(both_dir_shared);
    for (const auto& other_a : other_assemblies) {
      other_a->other_ref_coverage.push_back(both_dir_shared);
    }
  }
}

void calc_coverage::on_assembly(assembly_ptr a) {
  if (a->matches_reference) {
    return;
  }
  track_left_offset(a->left_offset);
  skip_to(a->left_offset - aoffset_t(m_options.seqset->max_read_len()));
  advance_to(a->left_offset);

  CHECK(m_ref_path_group);

  if (k_cov_dbg) {
    std::cout << "Adding variant at " << m_cur_offset << ": " << *a << "\n";
    m_ref_path_group->dump_debug_state();
  }

  auto sa = std::make_shared<shared_assembly>();
  sa->a = std::move(a);
  sa->c = this;

  for (auto matches : {&assembly::left_pair_matches, &assembly::right_pair_matches}) {
    if (k_cov_dbg) {
      if (!((*sa->a).*matches).empty()) {
        std::cout << "Variant looks for these read ids for pair coverage:";
      }
    }

    for (const auto& read_id : (*sa->a).*matches) {
      sa->pair_read_ids.insert(read_id);
      if (k_cov_dbg) {
        std::cout << " " << read_id;
      }
    }
    if (k_cov_dbg) {
      std::cout << "\n";
    }
  }

  auto var_dobj = make_unique<cov_tracker>();
  var_dobj->c = this;
  var_dobj->var_assembly = sa;
  size_t cov_size = sa->a->seq.size() + 1;
  size_t cov_elems_needed = coverage_accum::cov_elems_needed(cov_size);
  var_dobj->starts_and_ends.reset(new int[cov_elems_needed * 3]{});
  var_dobj->fwd_coverage.initialize(var_dobj->starts_and_ends.get() + 0 * cov_elems_needed,
                                    cov_size);
  var_dobj->rev_coverage.initialize(var_dobj->starts_and_ends.get() + 1 * cov_elems_needed,
                                    cov_size);
  var_dobj->pair_coverage.initialize(var_dobj->starts_and_ends.get() + 2 * cov_elems_needed,
                                     cov_size);
  for (const auto& active : m_active) {
    var_dobj->other_assemblies.push_back(active.second.first);
  }

  std::unique_ptr<path_group> var_pg = m_ref_path_group->split();
  var_pg->add_distant_object(std::move(var_dobj), sa->a->seq.size());
  if (k_cov_dbg) {
    std::cout << "After creating variant pg:\n";
    var_pg->dump_debug_state();
  }
  var_pg->add_sequence(sa->a->seq);

  if (k_cov_dbg) {
    std::cout << "After adding variant " << *sa->a << ", ref pg:\n";
    m_ref_path_group->dump_debug_state();
    std::cout << "Variant pg:\n";
    var_pg->dump_debug_state();
  }

  var_pg->flush();
  m_need_coverage_to = std::max<aoffset_t>(
      m_need_coverage_to, sa->a->right_offset + aoffset_t(m_options.seqset->max_read_len()));
  m_active.emplace(sa->a->right_offset, std::make_pair(sa, std::move(var_pg)));
  m_need_ref_dobj = true;
}

void calc_coverage::advance_to(aoffset_t offset) {
  //  std::cout << "Advancing to " << offset << "\n";
  while (m_cur_offset < offset) {
    advance_towards(offset);
    flush_sorted_to(m_cur_offset);
  }
}

void calc_coverage::add_ref_dobj(int ref_len) {
  if (k_cov_dbg) {
    std::cout << "Adding ref dobj, ref pg:\n";
    m_ref_path_group->dump_debug_state();
  }
  if (!m_active.empty()) {
    auto ref_dobj = make_unique<cov_tracker>();
    ref_dobj->c = this;
    size_t cov_size = ref_len + 1;
    size_t cov_elems_needed = coverage_accum::cov_elems_needed(cov_size);
    ref_dobj->starts_and_ends.reset(new int[cov_elems_needed * 2]{});
    ref_dobj->fwd_coverage.initialize(ref_dobj->starts_and_ends.get() + 0 * cov_elems_needed,
                                      cov_size);
    ref_dobj->rev_coverage.initialize(ref_dobj->starts_and_ends.get() + 1 * cov_elems_needed,
                                      cov_size);
    // TODO(nils): Get pair coverage for reference?
    for (const auto& active : m_active) {
      ref_dobj->other_assemblies.push_back(active.second.first);
    }
    m_ref_path_group->add_distant_object(std::move(ref_dobj), ref_len);
    if (k_cov_dbg) {
      std::cout << "After adding ref dobj, ref pg:\n";
      m_ref_path_group->dump_debug_state();
    }
  } else {
    if (k_cov_dbg) {
      std::cout << "Nothing active; skipping reference obj for this region\n";
    }
  }
  m_need_ref_dobj = false;
}

void calc_coverage::flush_active_to_here() {
  while (!m_active.empty() && m_active.begin()->first <= m_cur_offset) {
    if (m_need_ref_dobj) {
      // This may happen in the case of an insert; add coverage tracker before we rejoin.
      add_ref_dobj(0);
      CHECK(!m_need_ref_dobj);
    }
    CHECK(m_ref_path_group);
    if (k_cov_dbg) {
      std::cout << "Joining at " << m_cur_offset << ", ref pg:\n";
      m_ref_path_group->dump_debug_state();
      std::cout << "Joining at " << m_cur_offset << ", variant pg:\n";
      m_active.begin()->second.second->dump_debug_state();
    }
    m_ref_path_group->join(std::move(m_active.begin()->second.second));
    m_active.erase(m_active.begin());
    if (k_cov_dbg) {
      std::cout << "Post join at " << m_cur_offset << ", ref pg:\n";
      m_ref_path_group->dump_debug_state();
    }
  }
  CHECK(m_active.empty() || m_active.begin()->first > m_cur_offset);
}

void calc_coverage::skip_to(aoffset_t offset) {
  while (m_cur_offset < offset) {
    if (m_active.empty() && m_cur_offset >= m_need_coverage_to) {
      if (m_scaffold_it != m_options.scaffold->end() && m_scaffold_it.offset() < offset) {
        m_scaffold_it.skip_to(offset, "calc_coverage");
        m_cur_offset = m_scaffold_it.offset();
      } else {
        m_cur_offset = offset;
      }
      init_ref_pg();
    } else {
      //      std::cout << "Some still active; can't skip to " << offset << "\n";
      advance_towards(offset);
    }
  }
}

void calc_coverage::init_ref_pg() {
  if (k_cov_dbg && m_ref_path_group) {
    std::cout << "Resetting ref pg from:\n";
    m_ref_path_group->dump_debug_state();
  }
  m_ref_path_group = make_unique<path_group>(m_options.seqset->ctx_begin(),
                                             m_options.min_anchor_drop_overlap, this);
  m_ref_path_group->set_max_size(m_options.max_coverage_paths);
}

void calc_coverage::advance_towards(aoffset_t target_offset) {
  if (k_cov_dbg) {
    std::cout << "Advancing towards " << target_offset << " from " << m_cur_offset << "\n";
  }
  flush_active_to_here();

  CHECK_GT(target_offset, m_cur_offset);
  if (!m_active.empty() && target_offset > m_active.begin()->first) {
    target_offset = m_active.begin()->first;
  }

  CHECK_GT(target_offset, m_cur_offset);

  if (m_scaffold_it == m_options.scaffold->end()) {
    m_cur_offset = target_offset;
    return;
  } else if (m_scaffold_it.first_in_extent() && m_cur_offset < m_scaffold_it.offset()) {
    if (target_offset >= m_scaffold_it.offset()) {
      m_cur_offset = m_scaffold_it.offset();
      init_ref_pg();
      return;
    } else {
      m_cur_offset = target_offset;
      return;
    }
  }

  flush_active_to_here();

  if (m_scaffold_it.extent_end_offset() < target_offset) {
    target_offset = m_scaffold_it.extent_end_offset();
  }

  if (m_need_coverage_to > m_cur_offset && m_need_coverage_to < target_offset) {
    target_offset = m_need_coverage_to;
  }

  aoffset_t ref_len = target_offset - m_cur_offset;

  CHECK_GT(ref_len, 0);

  if (k_cov_dbg) {
    std::cout << "About to advance ref pg from " << m_cur_offset << " to " << target_offset
              << ", original ref pg is:\n";
    m_ref_path_group->dump_debug_state();
  }
  //  std::cout << "Making a ref obj of length " << ref_len << "\n";
  add_ref_dobj(ref_len);
  //  std::cout << "Adding bases from ref obj\n";
  if (k_cov_dbg) {
    std::cout << "ref pg is:\n";
    m_ref_path_group->dump_debug_state();
  }

  dna_sequence added_bases;
  while (m_cur_offset < target_offset) {
    CHECK(m_scaffold_it != m_options.scaffold->end());
    m_ref_path_group->add_base(*m_scaffold_it);
    if (k_cov_dbg) {
      added_bases.push_back(*m_scaffold_it);
    }
    //    //    std::cout << "Added " << *m_scaffold_it << "\n";
    //    m_ref_path_group->dump_debug_state();
    ++m_scaffold_it;
    ++m_cur_offset;
    if (m_cur_offset < target_offset) {
      CHECK_EQ(m_cur_offset, m_scaffold_it.offset());
      CHECK(m_scaffold_it != m_options.scaffold->end());
      CHECK(!m_scaffold_it.first_in_extent());
    }
  }
  m_ref_path_group->flush();

  if (k_cov_dbg) {
    std::cout << "Added " << added_bases << " to ref pg, which is now:\n";
    m_ref_path_group->dump_debug_state();
  }

  //  m_ref_path_group->dump_debug_state();
  //  std::cout << "Done making ref obj of length " << ref_len << "\n";
  //  m_ref_path_group->dump_debug_state();
}

void calc_coverage::on_seqset_entry(const seqset_range& r, path_group* pg) {
  if (r.size() < m_options.readmap->min_read_len()) {
    return;
  }
  auto reads = m_options.readmap->entry_to_index(r.begin());
  if (reads.first == reads.second) {
    return;
  }
  if (k_cov_dbg) {
    std::cout << "Found " << reads.second - reads.first << " reads for range "
              << r.sequence().rev_comp() << "\n";
  }
  dobj_visitor v(r, m_options.readmap, reads);
  pg->visit_distant_objects(r, v);
  if (k_cov_dbg) {
    std::cout << "Done applying coverage.\n";
  }
}

void calc_coverage::dobj_visitor::visit(path_group::distant_object* dobj, int distance) {
  cov_tracker* cov = dynamic_cast<cov_tracker*>(dobj);
  //  std::cout << "Visiting " << (void*)dobj << " with " << m_read_lens.size()
  //  << " read lengths, distance=" << distance << "\n";
  //  if (m_r.size() < distance) {
  //    return path_group::visit_result::CONSUME;
  //  }
  // if (distance > int(m_options.seqset->max_read_len())) {
  //   return path_group::visit_result::CONSUME;
  // }
  CHECK(cov);
  for (uint32_t read_id = m_reads.first; read_id != m_reads.second; ++read_id) {
    int read_len = m_readmap->get_readlength(read_id);
    if (read_len > int(m_r.size())) {
      continue;
    }
    bool is_forward = m_readmap->get_is_forward(read_id);

    auto& cov_accum = is_forward ? cov->fwd_coverage : cov->rev_coverage;

    aoffset_t read_start = aoffset_t(distance) - read_len + 1 + cov_accum.size() - 1;
    aoffset_t read_end = aoffset_t(distance) - 1 + cov_accum.size() - 1;

    if (k_cov_dbg) {
      std::cout << "Applying coverage start=" << read_start << " end=" << read_end
                << " (after distance=" << distance << ", cov size = " << cov_accum.size()
                << ") to ";
      if (cov->var_assembly) {
        std::cout << *cov->var_assembly->a << " with " << cov->other_assemblies.size()
                  << " other active";
      } else {
        std::cout << "a reference assembly with " << cov->other_assemblies.size() << " active";
      }
      std::cout << "\n";
    }

    if (cov_accum.add(read_start, read_end)) {
      if (cov->var_assembly) {
        if (cov->var_assembly->pair_read_ids.count(read_id)) {
          if (k_cov_dbg) {
            std::cout << "Pair matched; adding pair coverage\n";
          }
          CHECK(cov->pair_coverage.add(read_start, read_end));
        } else {
          if (k_cov_dbg) {
            std::cout << "Pair read_id=" << read_id << " did not match\n";
          }
        }
      }
    }
  }
}

calc_coverage::~calc_coverage() {
  skip_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
  m_ref_path_group->flush();
  m_ref_path_group.reset();
}

}  // namespace variants

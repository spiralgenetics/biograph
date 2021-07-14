#include "modules/variants/pop_tracer.h"

#include "modules/bio_base/readmap.h"
#include "modules/variants/assemble.h"
#include "modules/variants/ref_map.h"

namespace variants {

const char pop_tracer::k_pop_tracer_name[] = "POP";

static constexpr int k_trace_dbg = 0;

static std::set<uint32_t> g_debug_read_ids;
static std::set<seqset_range> g_debug_seqset_entries;

std::ostream& operator<<(std::ostream& os, const pop_tracer::entry& p) {
  os << "[" << p.left_offset << "," << p.right_offset << "): ~[" << p.start_limit << ","
     << p.end_limit << "] (orig len=";
  if (p.orig_r.valid()) {
    os << p.orig_r.size();
  } else {
    os << "invalid";
  }
  os << ") " << p.seq << " ";
  if (p.popped_r.valid()) {
    os << p.popped_r.sequence();
  } else {
    os << "invalid";
  }
  return os;
}

void pop_tracer::add_debug_read(uint32_t read_id) { g_debug_read_ids.insert(read_id); }

void pop_tracer::add_debug_seqset_entry(const seqset_range& r) { g_debug_seqset_entries.insert(r); }

void pop_tracer::clear_debug_reads() {
  g_debug_read_ids.clear();
  g_debug_seqset_entries.clear();
}

bool pop_tracer::range_needs_trace(const seqset_range& r) {
  for (const auto& debug_r : g_debug_seqset_entries) {
    if (debug_r.end() > r.begin() && debug_r.begin() < r.end()) {
      return true;
    }
  }
  return false;
}

bool pop_tracer::entry_needs_trace(const entry& e) {
  return range_needs_trace(e.orig_r) || range_needs_trace(e.popped_r) || e.trace_this;
}

pop_tracer::pop_tracer(const assemble_options& options) : m_options(options) {
  for (uint32_t read_id : g_debug_read_ids) {
    const readmap* rm = m_options.readmap;
    std::vector<uint32_t> expanded = {read_id, rm->get_rev_comp(read_id)};
    if (rm->has_mate(read_id) && false) {
      uint32_t rc = rm->get_rev_comp(read_id);
      expanded.push_back(rc);
      expanded.push_back(rm->get_rev_comp(rc));
    }
    for (uint32_t expand_read_id : expanded) {
      g_debug_seqset_entries.insert(rm->get_read_by_id(expand_read_id).get_seqset_entry());
    }
  }
}

pop_tracer::~pop_tracer() {}

void pop_tracer::add_read(uint32_t read_id, aoffset_t start_offset, aoffset_t limit_offset) {
  auto rd = m_options.readmap->get_read_by_id(read_id);
  auto r = rd.get_seqset_entry();

  std::shared_ptr<entry> new_entry = std::make_shared<entry>();
  new_entry->start_limit = start_offset;
  new_entry->end_limit = limit_offset;

  new_entry->orig_r = r;
  new_entry->popped_r = r;

  if (range_needs_trace(r) || k_trace_dbg > 3) {
    std::cout << "Provided range " << r.sequence() << " as a read, id " << read_id << " offset ["
              << start_offset << "," << limit_offset << ")\n";
  }

  add_entry_reads(*new_entry, r);
  m_fronts.emplace(r, new_entry);
}

void pop_tracer::add_reference_at(aoffset_t offset, const seqset_range& r) {
  if (r.size() < m_options.min_pop_overlap) {
    return;
  }

  bool trace_this = range_needs_trace(r);
  std::shared_ptr<entry> new_entry = std::make_shared<entry>();
  new_entry->left_offset = offset;
  new_entry->start_limit = offset;
  new_entry->end_limit = offset;
  new_entry->orig_r = r;
  new_entry->popped_r = r;
  new_entry->matches_reference = true;
  new_entry->trace_this = trace_this;
  add_entry_reads(*new_entry, r);
  m_fronts.emplace(r, new_entry);

  // if (!m_options.trace_ambiguous_ref) {
  //   if (m_options.rmap->has_potential_ambiguity(r)) {
  //     continue;
  //   }
  // }

  if (k_trace_dbg > 4 || trace_this) {
    std::cout << "Adding reference at " << offset << ": " << r.sequence().rev_comp()
              << "\nEntry: " << *new_entry << "\n";
  }

  //  m_poppers.emplace_back(std::move(new_entry));

  unsigned read_count = 0;

  for (const auto& rd : m_options.readmap->get_prefix_reads(r, m_options.min_overlap)) {
    if (!rd.has_mate()) {
      continue;
    }
    ++read_count;
    if (read_count > m_options.max_pairs_per_read) {
      if (k_trace_dbg > 4 || trace_this) {
        std::cout << "Too many reads; not checking any more\n";
      }
      break;
    }
    auto mate_rc = rd.get_mate().get_rev_comp();
    if (rd.size() < aoffset_t(m_options.min_pop_overlap)) {
      if (k_trace_dbg > 4 || trace_this) {
        std::cout << "Mate of read " << rd.get_read_id()
                  << " smaller than min_pop_overlap; ignoring\n";
      }
      continue;
    }

    aoffset_t start_offset = offset;
    aoffset_t limit_offset = offset;
    if (rd.is_original_orientation() == m_options.forward_pairs_face_inward) {
      // Faces inward.
      start_offset += rd.size() + m_options.min_pair_distance;
      limit_offset += rd.size() + m_options.max_pair_distance;
    } else {
      // Faces outwards
      start_offset -= rd.size() + m_options.max_pair_distance;
      limit_offset -= rd.size() + m_options.min_pair_distance;
    }

    add_read(mate_rc.get_read_id(), start_offset, limit_offset);
  }
}

void pop_tracer::add_reference(aoffset_t start, aoffset_t limit) {
  if (k_trace_dbg) {
    std::cout << "Adding reference with pairs [" << start << ", " << limit
              << ") to pop tracer, min_pop_overlap = " << m_options.min_pop_overlap << "\n";
  }
  // Reverse the scaffold, so the entries face right when we're doing
  // push_front_drop through it.  Then when we pop, we'll go towards
  // the more intuitive right.
  scaffold rc_scaffold = m_options.scaffold->rev_comp();
  aoffset_t scaffold_end = rc_scaffold.end_pos();
  scaffold::iterator next_it = rc_scaffold.begin();
  seqset_range r = m_options.seqset->ctx_begin();
  aoffset_t skip_pos = scaffold_end - limit - aoffset_t(m_options.readmap->max_read_len());
  if (skip_pos > next_it.offset()) {
    next_it.skip_to(skip_pos, "pop_tracer");
  }
  aoffset_t offset = scaffold_end - next_it.offset();
  while (next_it != rc_scaffold.end() && offset >= start) {
    if (next_it.first_in_extent()) {
      add_reference_at(offset, r);
      r = m_options.seqset->ctx_begin();
    }
    seqset_range next_r = r.push_front_drop((*next_it).complement());
    if (next_r.size() <= r.size()) {
      add_reference_at(offset, r);
    }

    r = next_r;
    if (next_it.first_in_extent()) {
      offset = scaffold_end - next_it.offset() - 1;
    } else {
      CHECK_EQ(offset, scaffold_end - next_it.offset());
      --offset;
    }

    ++next_it;
  }
  add_reference_at(offset, r);

  if (k_trace_dbg > 3) {
    std::cout << "Generated " << m_poppers.size() << " poppers and " << m_fronts.size()
              << " fronts\n";
    for (const auto& f : m_fronts) {
      std::cout << " F: " << *f.second << "\n";
    }
  }
}

void pop_tracer::add_anchor_drop(const assembly& a, bool right_anchor) {
  if (k_trace_dbg > 1 || assembly_needs_trace(a)) {
    std::cout << "IN: Pop tracer " << this << " adding anchor drop:" << a
              << " right_anchor=" << right_anchor << "\n";
  }

  aoffset_t offset = right_anchor ? a.right_offset : a.left_offset;

  // Add mates associated with these reads.
  for (uint32_t rc_read_id : a.rc_read_ids) {
    readmap::read rc_rd = m_options.readmap->get_read_by_id(rc_read_id);
    if (!rc_rd.has_mate()) {
      continue;
    }
    aoffset_t start_offset = offset - a.seq.size();
    aoffset_t limit_offset = offset + a.seq.size();
    if (right_anchor) {
      start_offset -= a.seq.size();
    } else {
      limit_offset += a.seq.size();
    }
    readmap::read mate = rc_rd.get_mate();
    // Mate faces in the forward direction.
    if (rc_rd.is_original_orientation() == m_options.forward_pairs_face_inward) {
      // rc_rd faces to the left, and this is inward towards the mate.
      start_offset -= rc_rd.size() + m_options.max_pair_distance;
    } else {
      limit_offset += rc_rd.size() + m_options.max_pair_distance;
    }
    add_read(mate.get_read_id(), start_offset, limit_offset);
  }

  // Add the half aligned section as an option to attach to.
  std::shared_ptr<entry> new_entry = std::make_shared<entry>();
  new_entry->start_limit =
      offset - aoffset_t(a.seq.size() + aoffset_t(m_options.max_pair_distance));
  new_entry->end_limit = offset + aoffset_t(a.seq.size() + m_options.max_pair_distance);

  if (right_anchor) {
    new_entry->right_offset = a.right_offset;
    new_entry->orig_r = m_options.seqset->ctx_begin();
    for (dna_base b : a.seq.rev_comp()) {
      new_entry->orig_r = new_entry->orig_r.push_front_drop(b.complement());
    }
    new_entry->popped_r = m_options.seqset->ctx_begin();
    new_entry->seq = a.seq;
    new_entry->trace_this = entry_needs_trace(*new_entry);
    if (k_trace_dbg > 1 || new_entry->trace_this) {
      std::cout << "Resultant entry from right half-aligned: " << *new_entry << "\n";
    }
  } else {
    new_entry->left_offset = a.left_offset;
    new_entry->orig_r = m_options.seqset->ctx_begin();
    for (dna_base b : a.seq.rev_comp()) {
      new_entry->orig_r = new_entry->orig_r.push_front_drop(b.complement());
    }
    new_entry->popped_r = m_options.seqset->ctx_begin();
    for (dna_base b : a.seq.rev_comp()) {
      seqset_range pushed = new_entry->popped_r.push_front(b.complement());
      if (!pushed.valid()) {
        break;
      }
      new_entry->popped_r = pushed;
    }
    CHECK_GE(a.seq.size(), new_entry->popped_r.size());
    new_entry->seq = a.seq.subseq(0, a.seq.size() - new_entry->popped_r.size());
    new_entry->trace_this = entry_needs_trace(*new_entry);
    if (k_trace_dbg > 1 || new_entry->trace_this) {
      std::cout << "Resultant entry from left half-aligned: " << *new_entry << "\n";
    }
    add_popper(new_entry);
  }
  seqset_range r = new_entry->orig_r;
  m_fronts.emplace(r, std::move(new_entry));
}

struct pop_tracer::popper_sorter {
  bool operator()(const std::shared_ptr<entry>& lhs, const std::shared_ptr<entry>& rhs) const {
    // Pop longer sizes first.
    return lhs->popped_r.size() < rhs->popped_r.size();
  }
};

void pop_tracer::assemble(assemble_pipeline_interface* output) {
  if (k_trace_dbg) {
    std::cout << "Starting pop trace assemble with " << m_fronts.size() << " fronts and "
              << m_poppers.size() << " poppers.  Poppers:\n";
    for (const auto& p : m_poppers) {
      std::cout << "  " << *p << "\n";
    }
  }

  std::make_heap(m_poppers.begin(), m_poppers.end(), popper_sorter());
  match_and_output_pass(output);
  while (!m_poppers.empty()) {
    pop_pass();
    match_and_output_pass(output);
  }
  if (k_trace_dbg) {
    std::cout << "Done pop trace assemble\n";
  }
}

struct pop_tracer::match_sorter {
  match_sorter() = delete;
  match_sorter(const match_sorter&) = default;
  match_sorter(const entry& p)
      : m_p(p), m_p_middle((m_p.start_limit + m_p.end_limit) / 2 + p.seq.size()) {}

  bool operator()(fronts_t::iterator lhs_it, fronts_t::iterator rhs_it) const {
    const entry& lhs = *lhs_it->second;
    const entry& rhs = *rhs_it->second;

    if (bool(lhs.left_offset) != bool(rhs.left_offset)) {
      // Reanchor to reference if possible!
      return bool(lhs.left_offset) > bool(rhs.left_offset);
    }

    if (lhs.orig_r.size() != rhs.orig_r.size()) {
      // More specific is better.
      return lhs.orig_r.size() > rhs.orig_r.size();
    }

    aoffset_t lhs_middle = (lhs.start_limit + lhs.end_limit) / 2;
    aoffset_t rhs_middle = (rhs.start_limit + rhs.end_limit) / 2;
    aoffset_t lhs_dist = abs(m_p_middle - lhs_middle);
    aoffset_t rhs_dist = abs(m_p_middle - rhs_middle);
    if (lhs_dist != rhs_dist) {
      return lhs_dist < rhs_dist;
    }

    if (lhs.seq.size() != rhs.seq.size()) {
      return rhs.seq.size() > rhs.seq.size();
    }

    aoffset_t lhs_span = lhs.end_limit - lhs.start_limit;
    aoffset_t rhs_span = rhs.end_limit - rhs.start_limit;
    if (lhs_span != rhs_span) {
      return lhs_span < rhs_span;
    }

    return false;
  }

  const entry& m_p;
  const aoffset_t m_p_middle;
};

void pop_tracer::pop_pass() {
  size_t orig_poppers_size = m_poppers.size();
  if (k_trace_dbg > 2) {
    std::cout << "Starting pop pass with " << m_poppers.size() << " poppers and " << m_fronts.size()
              << " fronts\n";
  }
  std::vector<std::shared_ptr<entry>> old_poppers;
  old_poppers.reserve(m_poppers.size());
  std::swap(old_poppers, m_poppers);

  while (!old_poppers.empty()) {
    std::pop_heap(old_poppers.begin(), old_poppers.end(), popper_sorter());
    std::shared_ptr<entry> p = std::move(old_poppers.back());
    old_poppers.pop_back();

    CHECK(!p->right_offset) << "Should not pop any more once we align the right side";
    CHECK_GE(p->popped_r.size(), m_options.min_pop_overlap);
    bool trace_this = entry_needs_trace(*p);
    seqset_range popped = p->popped_r.pop_front();
    trace_this = trace_this || range_needs_trace(popped);
    if (trace_this || k_trace_dbg > 4) {
      std::cout << "Popping popper: " << *p << "\n";
    }
    if (false) {
      if (popped.size() < m_options.min_pop_overlap) {
        if (popped.begin() + 1 == popped.end()) {
          // Trace a nonambiguous path
          uint64_t seqset_id = popped.begin();
          if (m_options.seqset->entry_size(seqset_id) > popped.size()) {
            popped = m_options.seqset->ctx_entry(seqset_id).truncate(popped.size() + 1);
          }
        }
      }
    }

    if (popped.size() < m_options.min_pop_overlap) {
      if (k_trace_dbg > 2 || trace_this) {
        std::cout << "Popper " << *p << " popped too much! Discarding\n";
      }
      continue;
    }

    p->seq.push_back(p->popped_r.front());
    p->popped_r = popped;
    add_entry_reads(*p, popped);
    p->end_limit++;
    add_popper(std::move(p));
  }
  if (k_trace_dbg) {
    if (m_poppers.size() != orig_poppers_size) {
      std::cout << "Pop pass decreased popper count from " << orig_poppers_size << " to " << m_poppers.size() << "\n";
    }
  }
}

void pop_tracer::add_popper(std::shared_ptr<entry> p) {
  if (p->left_offset) {
    CHECK_NE(m_options.scaffold->subscaffold_str(std::max<aoffset_t>(*p->left_offset, 0), 2), "NN") << *p << "\nRef region, -100:" << m_options.scaffold->subscaffold_str(*p->left_offset - 100, 200);
  }
  if (p->right_offset) {
    CHECK_NE(m_options.scaffold->subscaffold_str(std::max<aoffset_t>(*p->right_offset, 0), 2), "NN") << *p;
  }
  m_poppers.emplace_back(std::move(p));
  std::push_heap(m_poppers.begin(), m_poppers.end(), popper_sorter());
}

void pop_tracer::add_entry_reads(entry& p, const seqset_range& r) {
  size_t count = 0;
  for (const auto& rd : m_options.readmap->get_prefix_reads(r, m_options.min_overlap)) {
    ++count;
    if (count > m_options.max_pairs_per_read) {
      break;
    }
    p.seen_read_ids.push_back(rd.get_rev_comp().get_read_id());
  }
}

void pop_tracer::match_and_output_pass(assemble_pipeline_interface* output) {
  size_t orig_poppers_size = m_poppers.size();
  if (k_trace_dbg > 2) {
    std::cout << "Starting match and output pass with " << m_poppers.size() << " poppers and "
              << m_fronts.size() << " fronts\n";
  }
  std::vector<std::shared_ptr<entry>> old_poppers;
  old_poppers.reserve(m_poppers.size());
  std::swap(old_poppers, m_poppers);
  for (auto& p : old_poppers) {
    bool trace_this = entry_needs_trace(*p);
    if (trace_this || k_trace_dbg > 4) {
      std::cout << "Matching popper: " << *p << "\n";
    }

    bool found_match = false;

    std::vector<fronts_t::iterator> matches;

    // See if we can find a match for this!
    for (auto it = m_fronts.lower_bound(p->popped_r);
         it != m_fronts.end() && it->second->orig_r.begin() < p->popped_r.end(); ++it) {
      const entry& front = *it->second;
      bool trace_this_match = trace_this;

      if (trace_this || entry_needs_trace(front)) {
        trace_this_match = true;
        std::cout << "Does Popper " << *p << " match " << front << "?\n";
      }

      if (&front == p.get()) {
        if (trace_this_match) {
          std::cout << "Match failed: loop\n";
        }
        // no loops!
        continue;
      }

      if (front.orig_r.end() > p->popped_r.end()) {
        // Not a prefix.
        if (trace_this_match) {
          std::cout << "Not actually a prefix\n";
        }
        continue;
      }

      if (p->end_limit + aoffset_t(p->seq.size() + m_options.pop_tracer_offset_slop) <
          front.start_limit) {
        // Not close enough to use.
        if (trace_this_match) {
          std::cout << "Out of range (1)\n";
        }
        continue;
      }

      if (front.end_limit + aoffset_t(p->seq.size() + m_options.pop_tracer_offset_slop) <
          p->start_limit) {
        // Not close enough to use.
        if (trace_this_match) {
          std::cout << "Out of range (2)\n";
        }
        continue;
      }

      if (front.right_offset && p->left_offset && *front.right_offset < *p->left_offset) {
        if (trace_this_match) {
          std::cout << "Misordered!?\n";
        }
        continue;
      }
      if (front.left_offset && p->left_offset && *front.left_offset < *p->left_offset) {
        if (trace_this_match) {
          std::cout << "Misordered left offset?\n";
        }
        continue;
      }

      if (trace_this_match) {
        std::cout << "Adding to match list! Enabling tracing for this whole popper.\n";
        trace_this = true;
      }

      if (k_trace_dbg > 3 || trace_this) {
        std::cout << "Popper " << *p << " #" << matches.size() << " match with front: " << front
                  << "\n";
      }
      matches.push_back(it);
    }

    if ((k_trace_dbg > 2 || trace_this) && matches.size() > 1) {
      std::cout << "Popper " << *p << " ambiguously found " << matches.size() << " matches\n";
    }

    std::sort(matches.begin(), matches.end(), match_sorter(*p));

    for (auto it : matches) {
      const entry& front = *it->second;

      found_match = true;

      p->trace_this = p->trace_this || front.trace_this;

      // Looks like we got a matching front!
      if (front.left_offset) {
        // Front is aligned to reference.
        p->right_offset = front.left_offset;
        if (k_trace_dbg > 3 || trace_this) {
          std::cout << "Popper " << *p << " joined to reference on right at " << *front.left_offset
                    << ", not adding back to queue\n";
        }

        if (p->left_offset) {
          // Successfully attached on both sides!  Output assembly if this isn't completely
          // reference-only.
          if (k_trace_dbg > 3 || trace_this) {
            std::cout << "Outputting assembly for " << *p << "\n";
          }
          output_assembly(output, *p);
        } else {
          if (k_trace_dbg > 2 || trace_this) {
            std::cout << "Popper " << *p << " found right alignment; waiting for left.\n";
          }
        }
      } else {
        // Merge in.
        p->start_limit = std::max(p->start_limit, front.start_limit - aoffset_t(p->seq.size()));
        p->end_limit = std::min(p->end_limit, front.end_limit - aoffset_t(p->seq.size()));
        p->popped_r = front.popped_r;
        p->seq += front.seq;
        p->trace_this = p->trace_this || front.trace_this;
        p->seen_read_ids.insert(p->seen_read_ids.end(), front.seen_read_ids.begin(), front.seen_read_ids.end());
        if (k_trace_dbg > 2 || trace_this) {
          std::cout << "Popper merged with " << front << ", resulting in:\n" << *p << "\n";
        }
        if (front.right_offset) {
          p->right_offset = front.right_offset;
          if (k_trace_dbg > 2 || trace_this) {
            std::cout << "Popper " << *p << " joined right alignment; still waiting for left\n";
          }
          if (p->left_offset) {
            // Alignment successful!
            if (k_trace_dbg > 2 || trace_this) {
              std::cout << "Popper " << *p << " generating assembly\n";
            }
            output_assembly(output, *p);
          }
        } else {
          if (trace_this) {
            std::cout << "Popper " << *p << " still needs popping after merging\n";
          }
          // Re-add the merged entry to keep searching for something to align to on the right.
          m_poppers.emplace_back(std::move(p));
        }
      }

      if (front.left_offset) {
        // Leave reference sections so anyone can align to them.
      } else {
        // Otherwise, consume this right side in order to save space.
        // This also makes it so we can't loop infinitely.
        m_fronts.erase(it);
      }
      break;
    }
    if (!found_match) {
      // Try again after popping more.
      m_poppers.emplace_back(std::move(p));
    }
  }

  if (k_trace_dbg) {
    if (m_poppers.size() != orig_poppers_size) {
      std::cout << "Match and output pass decreased popper count from " << orig_poppers_size << " to " << m_poppers.size() << ":\n";
      for (const auto& p : m_poppers) {
        std::cout << "  " << *p << "\n";
      }
    }
  }
}

void pop_tracer::output_assembly(assemble_pipeline_interface* output, const entry& p) {
  CHECK(p.left_offset);
  CHECK(p.right_offset);

  assembly_ptr a = make_unique<assembly>();
  a->tags.insert(k_pop_tracer_name);
  a->assembly_id = allocate_assembly_id();
  a->left_offset = *p.left_offset;
  a->right_offset = *p.right_offset;
  a->seq = p.seq;
  a->rc_read_ids.insert(p.seen_read_ids.begin(), p.seen_read_ids.end());

  if (a->left_offset == a->right_offset && a->seq.size() == 0) {
    if (k_trace_dbg) {
      std::cout << "Null assembly generated from " << p << "\n";
    }
    return;
  }

  if (a->seq.size() <= p.orig_r.size() &&
      aoffset_t(p.seq.size()) == a->right_offset - a->left_offset) {
    if (k_trace_dbg) {
      a->matches_reference = true;
      std::cout << "Reference assembly: " << *a << "\n";
    }
  }
  if (k_trace_dbg > 4) {
    std::cout << "read ids:\n";
    for (uint32_t read_id : a->rc_read_ids) {
      std::cout << " " << read_id << " "
                << m_options.readmap->get_read_by_id(read_id).get_seqset_entry().sequence() << "\n";
      ;
    }
  }
  if (assembly_needs_trace(*a)) {
    std::cout << "OUT: pop_tracer " << this << " produced " << (void *)a.get() << ": " << *a << "\n";
  }
  output->add(std::move(a));
}

}  // namespace variants

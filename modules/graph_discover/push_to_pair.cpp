#include "modules/graph_discover/push_to_pair.h"

#include <boost/range/adaptor/reversed.hpp>
#include <chrono>
#include <ctime>
#include <random>

#include "modules/bio_base/readmap.h"
#include "modules/io/parallel.h"

namespace variants {

constexpr bool k_dbg = false;
constexpr bool k_check_reads = false;

std::string push_to_pair_discover::push_active_assembly::to_string() const {
  std::stringstream os;
  os << active_assembly::to_string() << ", " << mates.size() << " mates =";
  for (const auto& mate : mates) {
    os << " " << mate.sequence();
  }
  return os.str();
}

void push_to_pair_discover::flush() {
  graph_discover::flush();
  CHECK(m_rc_entry_anchors.empty());
}

push_to_pair_discover::push_to_pair_discover(const assemble_options& options,
                                             const std::string& tag, pipeline_step_t output)
    : graph_discover(options, std::move(output)), m_tag(tag) {
  CHECK(opts().readmap);
}

push_to_pair_discover::~push_to_pair_discover() {
  CHECK(m_mates.empty());
  CHECK(m_mate_expiry.empty());
}

void push_to_pair_discover::on_walk(active_assembly* generic_act) {
  push_active_assembly* act = dynamic_cast<push_active_assembly*>(generic_act);
  CHECK(act);

  act->mates = act->a->rc_seqset_entries.mates();
  aoffset_t len = act->a->seq.size();

  if (act->a->right_offset) {
    for (const auto& e : act->a->rc_seqset_entries.entries()) {
      aoffset_t offset = e.first;
      if (offset == len) {
        continue;
      }
      for (const auto& r : e.second) {
        if (r.size() < opts().min_overlap) {
          continue;
        }
        anchor_info info;
        info.truncated = r.truncate(opts().min_overlap);
        info.orig = r;
        info.anchor.act = act;
        info.anchor.offset = len - e.first;
        act->rc_anchors.emplace_back(std::move(info));
      }
    }
  }

  if (k_dbg) {
    std::cerr << "Found " << act->mates.size() << " mates for " << *act->a << ":\n";
    for (const auto& mate : act->mates) {
      std::cerr << " " << mate.sequence() << "\n";
    }
  }
}

void push_to_pair_discover::on_advance_trace(aoffset_t offset) {
  auto first_exp = m_mate_expiry.begin();
  while (first_exp != m_mate_expiry.end() && first_exp->first < offset) {
    if (k_dbg) {
      std::cerr << "Removing mates at " << first_exp->first << "\n";
    }
    for (const auto& mate : first_exp->second) {
      auto it = m_mates.find(mate);
      CHECK(it != m_mates.end());
      CHECK_GT(it->second, 0);
      --it->second;
      if (!it->second) {
        m_mates.erase(it);
      }
    }
    first_exp = m_mate_expiry.erase(first_exp);
  }
  if (k_dbg) {
    if (first_exp != m_mate_expiry.end()) {
      std::cerr << "Next mate expiry at " << first_exp->first << "\n";
    } else {
      std::cerr << "No mate expiries pending\n";
    }
  }
}

void push_to_pair_discover::on_readahead(const active_assembly* generic_act) {
  const push_active_assembly* act = dynamic_cast<const push_active_assembly*>(generic_act);
  CHECK(act);

  aoffset_t right_asm_offset = max(act->a->left_offset, act->a->right_offset);
  aoffset_t exp_offset = right_asm_offset + opts().max_pair_distance;
  auto& expiry = m_mate_expiry[exp_offset];

  if (k_dbg) {
    std::cerr << "Adding " << act->mates.size() << " mates from " << *act->a << " to expire at "
              << exp_offset << "(" << right_asm_offset << " + " << opts().max_pair_distance
              << ")\n";
  }

  for (const auto& act_info : act->rc_anchors) {
    auto& info = m_rc_entry_anchors[act_info.truncated];
    info.rc_starts.emplace(act_info.orig, act_info.anchor);
  }

  for (const auto& mate : act->mates) {
    if (expiry.insert(mate).second) {
      ++m_mates[mate];
    }
  }
}

void push_to_pair_discover::on_readahead_done(const active_assembly* generic_act) {
  const push_active_assembly* act = dynamic_cast<const push_active_assembly*>(generic_act);
  CHECK(act);

  if (act->a->right_offset) {
    for (const auto& act_info : act->rc_anchors) {
      auto info_it = m_rc_entry_anchors.find(act_info.truncated);
      CHECK(info_it != m_rc_entry_anchors.end());
      auto& info = info_it->second;

      auto anchors = info.rc_starts.equal_range(act_info.orig);
      for (auto it = anchors.first;; ++it) {
        CHECK(it != anchors.second) << "Could not find entry anchor to delete";
        if (it->first != act_info.orig) {
          continue;
        }
        const auto& p = it->second;
        if (p.act != act_info.anchor.act || p.offset != act_info.anchor.offset) {
          continue;
        }
        info.rc_starts.erase(it);
        break;
      }

      if (info.rc_starts.empty()) {
        m_rc_entry_anchors.erase(info_it);
      }
    }
  }
}

seqset_range_set push_to_pair_discover::trace_one_base(const seqset_range_set& rs,
                                                       dna_sequence* seq, int* bases_since_read) {
  dna_base_array<seqset_range_set> next_rs;
  dna_base_array<seqset_range_set> next_rs_drop;
  for (dna_base rc_b : dna_bases()) {
    for (const auto& r : rs) {
      seqset_range next_r = r.push_front_drop(rc_b, opts().min_overlap + *bases_since_read);

      if (!next_r.valid()) {
        continue;
      }

      bool drop = next_r.size() != (r.size() + 1);

      if (drop) {
        next_rs_drop[rc_b].insert(next_r);
      } else {
        next_rs[rc_b].insert(next_r);
      }
    }
  }

  unsigned best_len = 0;
  dna_base best_rc_b;

  for (dna_base rc_b : dna_bases()) {
    unsigned next_len = 0;
    for (const auto& r : next_rs[rc_b]) {
      if (r.size() > next_len) {
        next_len = r.size();
      }
    }
    for (const auto& r : next_rs_drop[rc_b]) {
      if (r.size() > next_len) {
        next_len = r.size();
      }
    }
    if (next_len > best_len) {
      best_len = next_len;
      best_rc_b = rc_b;
    }
  }
  if (!best_len) {
    return {};
  }

  auto& best_rs = next_rs[best_rc_b];
  auto& best_rs_drop = next_rs_drop[best_rc_b];
  dna_base best_b = best_rc_b.complement();

  bool found_read = false;
  if (!k_check_reads) {
    found_read = true;
  }
  for (const auto& r : best_rs) {
    if (!found_read) {
      if (opts().readmap->get_longest_prefix_read_id(r)) {
        found_read = true;
      }
    }
  }
  for (const auto& r : best_rs_drop) {
    if (!m_seen_ranges.insert(r).second) {
      // best is a duplicate.  stop tracing here.
      return {};
    }

    if (!found_read) {
      if (opts().readmap->get_longest_prefix_read_id(r)) {
        found_read = true;
      }
    }
  }
  if (found_read) {
    *bases_since_read = 0;
  }
  (*seq).push_back(best_b);

  best_rs.insert(best_rs_drop.begin(), best_rs_drop.end());
  seqset_set_dedup_prefixes(best_rs);
  return std::move(best_rs);
}

void push_to_pair_discover::on_trace(const graph_discover::active_assembly* generic_act) {
  const push_active_assembly* act = dynamic_cast<const push_active_assembly*>(generic_act);
  CHECK(act);

  if (act->a->right_offset) {
    // This assembly is already anchored; don't extend it.
    return;
  }

  aoffset_t min_overlap = opts().min_overlap;
  if (k_dbg) {
    std::cerr << "Push to pair assembly " << *act->a << ", min overlap = " << min_overlap << "\n";
  }

  // Start at the right end of the assembly
  seqset_range_set rs = act->a->rc_seqset_entries.starts();

  // Make sure we don't loop back on ourself.
  m_seen_ranges.clear();

  dna_sequence seq;
  if (k_dbg) {
    std::cerr << "Searching from " << rs << "\n";
  }

  absl::btree_map<aoffset_t, seqset_range_set> path_entries;
  for (auto elem : act->a->rc_seqset_entries.entries()) {
    path_entries[act->a->seq.size() - elem.first] = std::move(elem.second);
  }

  aoffset_t best_abs_svlen = std::numeric_limits<aoffset_t>::max();
  assembly_ptr best_abs_svlen_a;
  aoffset_t best_shared_bases = 0;
  assembly_ptr best_shared_bases_a;
  int bases_since_read = 0;
  while (!rs.empty()) {
    seqset_range_set next_rs = trace_one_base(rs, &seq, &bases_since_read);

    if (next_rs.empty()) {
      // Nothing found
      if (k_dbg) {
        std::cerr << "Nothing found pushing here\n";
      }
      path_entries[act->a->seq.size() + seq.size()] = rs;
      break;
    }

    if (k_dbg) {
      std::cerr << "Found next, seen size = " << m_seen_ranges.size() << ":\n";
      for (const auto& r : next_rs) {
        std::cerr << " " << r.sequence() << "len=" << r.size() << "\n";
      }
    }

    rs = std::move(next_rs);

    bool found_anchor = false;
    bool found_mate = false;
    bool found_full_anchor = false;
    for (const auto& r : rs) {
      if (m_mates.count(r)) {
        if (k_dbg) {
          std::cerr << " " << r.sequence() << " has mate in other assembly\n";
        }
        found_mate = true;
      } else {
        if (k_dbg) {
          std::cerr << " " << r.sequence() << " has no mate\n";
        }
      }
    }

    aoffset_t max_size = 0;
    for (const auto& r : rs) {
      max_size = std::max<aoffset_t>(r.size(), max_size);
    }
    aoffset_t shared_bases_here =
        save_anchors(act, rs, seq, path_entries, best_abs_svlen, best_abs_svlen_a,
                     best_shared_bases, best_shared_bases_a);
    if (shared_bases_here) {
      if (k_dbg) {
        std::cerr << " " << rs << " saved anchors\n";
      }
      found_anchor = true;
      if (shared_bases_here == max_size) {
        found_full_anchor = true;
      }
    }

    if (found_mate || found_anchor) {
      if (k_dbg) {
        std::cerr << "Found mate or anchor!  Emitting partial.\n";
      }
      path_entries[act->a->seq.size() + seq.size()] = rs;
    }

    if (found_full_anchor) {
      // Emitted an anchor; don't trace anymore.
      break;
    }
  }

  if (best_shared_bases_a) {
    sort_and_output(std::move(best_shared_bases_a));
  }
  if (best_abs_svlen_a) {
    sort_and_output(std::move(best_abs_svlen_a));
  }

  if (!path_entries.empty()) {
    sort_and_output(discover_extend_right(
        act, act->a->seq.size(), seq, m_tag,
        path_entries_to_rc_path(path_entries, act->a->seq.size() + seq.size())));
  }
}

class push_to_pair_discover::entry_anchor_adder {
 public:
  entry_anchor_adder(std::vector<potential_entry_anchor>& anchors, aoffset_t approx_ref_offset)
      : m_anchors(anchors), m_approx_ref_offset(approx_ref_offset) {}

  void add(seqset_range orig_r, seqset_range r, const potential_anchor& anchor) {
    aoffset_t anchor_offset =
        anchor.act->a->right_offset - anchor.act->a->seq.size() + anchor.offset;
    aoffset_t abs_svlen = abs(anchor_offset - m_approx_ref_offset);
    if (abs_svlen > m_best_abs_svlen) {
      return;
    }
    m_best_abs_svlen = abs_svlen;

    potential_entry_anchor eanchor;
    (potential_anchor&)eanchor = anchor;
    eanchor.r = r;
    eanchor.shared_bases = eanchor.r.shared_prefix_length(orig_r);
    m_anchors.push_back(eanchor);
  }

 public:
  std::vector<potential_entry_anchor>& m_anchors;
  aoffset_t m_approx_ref_offset;
  aoffset_t m_best_abs_svlen = std::numeric_limits<aoffset_t>::max();
};

struct push_to_pair_discover::entry_anchor_search {
  rc_starts_t::const_iterator begin;
  rc_starts_t::const_iterator rev_it;
  rc_starts_t::const_iterator fwd_it;
  rc_starts_t::const_iterator end;
  seqset_range r;
};

std::vector<push_to_pair_discover::potential_entry_anchor> push_to_pair_discover::get_entry_anchors(
    const seqset_range_set& rs, aoffset_t min_overlap, aoffset_t approx_ref_offset) const {
  std::vector<potential_entry_anchor> anchors;
  CHECK_GT(min_overlap, 1);

  std::vector<entry_anchor_search> searches;

  for (const auto& r : rs) {
    if (aoffset_t(r.size()) < min_overlap) {
      continue;
    }
    auto it = m_rc_entry_anchors.upper_bound(r);
    if (it != m_rc_entry_anchors.begin()) {
      --it;

      if (it->first.end() <= r.begin()) {
        // None found
        continue;
      }

      const auto& info = it->second;

      entry_anchor_search s;
      s.begin = info.rc_starts.begin();
      s.end = info.rc_starts.end();
      s.fwd_it = info.rc_starts.upper_bound(r);
      s.rev_it = s.fwd_it;
      s.r = r;

      searches.push_back(s);
    }
  }

  entry_anchor_adder adder(anchors, approx_ref_offset);
  while (!searches.empty()) {
    std::vector<entry_anchor_search> new_searches;
    for (auto& s : searches) {
      bool more = false;
      if (s.rev_it != s.begin) {
        --s.rev_it;
        adder.add(s.r, s.rev_it->first, s.rev_it->second);
        more = true;
      }
      if (s.fwd_it != s.end) {
        adder.add(s.r, s.fwd_it->first, s.fwd_it->second);
        ++s.fwd_it;
        more = true;
      }

      if (more) {
        new_searches.push_back(s);
      }
    }
    searches = std::move(new_searches);
  }

  return anchors;
}

seqset_path push_to_pair_discover::path_entries_to_rc_path(
    const absl::btree_map<aoffset_t, seqset_range_set>& path_entries, aoffset_t end_pos) {
  CHECK(!path_entries.empty());
  CHECK_EQ(path_entries.begin()->first, 0);
  CHECK_LE(path_entries.rbegin()->first, end_pos);
  seqset_path rc_path;
  for (const auto& e : path_entries) {
    rc_path.add(end_pos - e.first, e.second);
  }
  return rc_path;
}

aoffset_t push_to_pair_discover::save_anchors(
    const active_assembly* act, const seqset_range_set& rs, dna_slice seq,
    const absl::btree_map<aoffset_t, seqset_range_set>& path_entries, aoffset_t& best_abs_svlen,
    assembly_ptr& best_abs_svlen_a, aoffset_t& best_shared_bases,
    assembly_ptr& best_shared_bases_a) {
  CHECK(!path_entries.empty());
  CHECK_LE(path_entries.rbegin()->first, act->a->seq.size() + seq.size());
  aoffset_t ref_offset = act->a->left_offset + act->a->seq.size() + seq.size();
  auto anchors = get_entry_anchors(rs, opts().min_overlap, ref_offset);
  if (anchors.empty()) {
    return 0;
  }
  auto abs_svlen = [ref_offset](const potential_entry_anchor& anchor) {
    aoffset_t anchor_offset =
        anchor.act->a->right_offset - anchor.act->a->seq.size() + anchor.offset;
    return abs(anchor_offset - ref_offset);
  };

  if (k_dbg) {
    std::cerr << anchors.size() << " potential anchors found\n";
  }

  seqset_path new_rc_path = path_entries_to_rc_path(path_entries, act->a->seq.size() + seq.size());
  new_rc_path.add(act->a->seq.size() + seq.size(), rs);

  std::sort(anchors.begin(), anchors.end(),
            [abs_svlen](const potential_entry_anchor& a, const potential_entry_anchor& b) {
              if (a.shared_bases != b.shared_bases) {
                return a.shared_bases > b.shared_bases;
              }
              auto a_svlen = abs_svlen(a);
              auto b_svlen = abs_svlen(b);
              if (a_svlen != b_svlen) {
                return a_svlen < b_svlen;
              }
              return false;
            });

  aoffset_t shared_bases_here = 0;
  for (const auto& anchor : anchors) {
    if (!anchor.act->a->right_offset) {
      continue;
    }
    aoffset_t cur_abs_svlen = abs_svlen(anchor);

    if (cur_abs_svlen < best_abs_svlen || anchor.shared_bases > best_shared_bases) {
      if (cur_abs_svlen < best_abs_svlen) {
        best_abs_svlen = cur_abs_svlen;
        best_abs_svlen_a = discover_anchor(act, act->a->seq.size() /* extending from the end */,
                                           seq, anchor, m_tag, new_rc_path);
      }
      if (anchor.shared_bases > best_shared_bases) {
        if (!shared_bases_here) {
          shared_bases_here = anchor.shared_bases;
        }
        best_shared_bases = anchor.shared_bases;
        best_shared_bases_a = discover_anchor(act, act->a->seq.size() /* extending from the end */,
                                              seq, anchor, m_tag, new_rc_path);
      }
    }
  }
  return shared_bases_here;
}

}  // namespace variants

#include "modules/variants/pair_cov.h"

#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "absl/container/btree_map.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"

namespace variants {

namespace {
constexpr int k_cov_debug = 0;
constexpr bool k_extended_stats = false;
constexpr bool k_report_progress_enabled = false;
constexpr int k_report_seconds = 300;
}  // namespace

struct pair_cov::result : private boost::noncopyable, public std::enable_shared_from_this<result> {
  assembly_ptr a;

  // Same as a->right_offset and a->left_offset; copied here to avoid
  // a pointer dereference since we look at these often.
  aoffset_t left_offset;
  aoffset_t right_offset;

  // seq.size() bases to the left of right_offset.
  aoffset_t right_anchored_left_offset;

  // Coverage entries that have been found to have pairing support.
  read_coverage_set pair_supported_reads;

  // Current sv adjustment bounds.
  // E.g. closed_hull(sv_adjusts.end().size_changes) or [0,1).
  interval_set_t cur_sv_adjust = closed_interval_set(0);
};

struct pair_cov::pair_cov_pg {
  // Pair table containing information from the assembly that's been processed in this path.
  pair_table pending_pair_table;

  result_ptr r;

  absl::flat_hash_map<multi_mid,
                      boost::container::flat_map<aoffset_t /* read offset */, read_id_mask_t>,
                      multi_mid_hash>
      pending_saved_reads;
};

std::string pair_cov::dump_read_id_mask(read_id_mask_t mask) const {
  std::stringstream out;

  for (unsigned i = 0; i != sizeof(read_id_mask_t) * 8; ++i) {
    out << (mask & 1);
    mask >>= 1;
  }
  return out.str();
}

std::string pair_cov::dump_pair_table(const pair_table& pt) const {
  std::stringstream out;

  if (k_cov_debug < 2) {
    return "";
  }
  for (const auto& rd_and_pte : pt.entries) {
    const multi_mid& mm = rd_and_pte.first;
    const auto& pte = rd_and_pte.second;
    out << "    " << dump_multi_mid(mm) << ", " << pte.results.size() << " read masks:\n";
    for (const auto& read_ids_and_roffs : pte.results) {
      const auto& masks = read_ids_and_roffs.first;
      const auto& roffs = read_ids_and_roffs.second;
      out << "      total=" << dump_read_id_mask(masks.total)
          << " pending=" << dump_read_id_mask(masks.pending) << ", " << roffs.size()
          << " results:\n";
      for (const auto& roff : roffs) {
        std::cout << "        read_start=" << roff.read_start << " in " << *roff.r->a << "\n";
      }
    }
  }
  return out.str();
}

std::string pair_cov::dump(pair_cov_pg* pg) const {
  std::stringstream out;

  out << "PairCov, main " << m_main_pair_table.entries.size() << " entries:\n";
  out << dump_pair_table(m_main_pair_table);
  out << "Pending, " << pg->pending_pair_table.entries.size() << " entries:\n";
  out << dump_pair_table(pg->pending_pair_table);
  return out.str();
}

pair_cov::pair_cov(const assemble_options& opts, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_opts(opts) {}

pair_cov::~pair_cov() {
  if (k_extended_stats) {
    std::cout << "Pair cov finishing up\n";
  }
  advance_ref_to(std::numeric_limits<aoffset_t>::max());
  flush_old();
  if (k_extended_stats) {
    std::cout << "Pair cov finished\n";
  }
  CHECK(m_active.empty());
  CHECK(m_cur_inserts.empty());
  CHECK(m_cur_non_inserts.empty());
  CHECK(m_pending_results.empty());
  CHECK(m_main_pair_table.entries.empty());
}

void pair_cov::advance_ref_to(aoffset_t offset) {
  while (m_cur_ref_offset < offset) {
    advance_ref_towards(offset);
    flush_sorted_to(m_cur_ref_offset);
  }
}

void pair_cov::advance_ref_towards(aoffset_t target_offset) {
  flush_active_to_here();

  if (!m_active.empty() && target_offset > m_active.begin()->first) {
    target_offset = m_active.begin()->first;
  }

  CHECK_LT(m_cur_ref_offset, target_offset);
  CHECK(m_cur_inserts.empty());
  CHECK(m_cur_non_inserts.empty());
  m_cur_ref_offset = target_offset;
}

void pair_cov::report_progress() {
  if (!k_report_progress_enabled) {
    return;
  }

  if (m_last_report == 0) {
    m_last_report = time(0);
    m_last_report_offset = m_cur_ref_offset;
    m_last_report_asms = m_cur_asm_count;
    return;
  }

  time_t now = time(0);
  if (m_last_report == 0) {
    m_last_report = now - (k_report_seconds / 2);
  } else if (m_last_report + k_report_seconds < now) {
    int elapsed = now - m_last_report;
    aoffset_t bases = m_cur_ref_offset - m_last_report_offset;
    int asms = m_cur_asm_count - m_last_report_asms;
    std::cout << printstring("%.2f", asms * 1. / elapsed) << " asm/sec, "
              << printstring("%.2f", bases * 1. / elapsed)
              << " bases/sec, offset=" << m_cur_ref_offset << " elapsed=" << elapsed
              << " bases=" << bases << " cur pend=" << m_pending_results.size()
              << " active=" << m_active.size() << "\n";

    m_last_report = now;
    m_last_report_offset = m_cur_ref_offset;
    m_last_report_asms = m_cur_asm_count;
  }
}

void pair_cov::flush_active_to_here() {
  if (k_cov_debug) {
    std::cout << "Flushing active to " << m_cur_ref_offset << "\n";
  }

  std::vector<std::unique_ptr<pair_cov_pg>> rejoins;
  while (!m_active.empty() && m_active.begin()->first == m_cur_ref_offset) {
    if (k_cov_debug) {
      std::cout << "Joining from assembly ending at " << m_cur_ref_offset << ": "
                << dump(m_active.begin()->second.get()) << "\n";
    }
    std::unique_ptr<pair_cov_pg> pg = std::move(m_active.begin()->second);
    rejoins.emplace_back(std::move(pg));
    m_active.erase(m_active.begin());
  }
  if (!m_active.empty()) {
    CHECK_GT(m_active.begin()->first, m_cur_ref_offset);
  }

  if (!rejoins.empty()) {
    join(std::move(rejoins));
  }

  update_adjusts();

  if (!m_cur_inserts.empty()) {
    // Split out a path group for each insert.
    std::vector<std::unique_ptr<pair_cov_pg>> insert_rejoins;
    for (auto& a : m_cur_inserts) {
      std::unique_ptr<pair_cov_pg> pg = make_unique<pair_cov_pg>();
      add_assembly(pg.get(), std::move(a));
      insert_rejoins.emplace_back(std::move(pg));
    }
    m_cur_inserts.clear();
    join(std::move(insert_rejoins));
  }

  update_adjusts();

  if (!m_cur_non_inserts.empty()) {
    for (auto& a : m_cur_non_inserts) {
      aoffset_t right_offset = a->right_offset;

      std::unique_ptr<pair_cov_pg> pg = make_unique<pair_cov_pg>();
      add_assembly(pg.get(), std::move(a));
      m_active.emplace(right_offset, std::move(pg));
    }
    m_cur_non_inserts.clear();
  }

  if (k_cov_debug) {
    std::cout << "Done flushing active to " << m_cur_ref_offset << "\n";
  }

  if (m_cur_ref_offset > m_next_flush_old) {
    m_next_flush_old = m_cur_ref_offset + aoffset_t(m_opts.max_pair_distance);
    flush_old();
    report_progress();
  }
}

void pair_cov::on_assembly(assembly_ptr a) {
  if (k_cov_debug) {
    std::cout << "Pair_cov received assembly: " << *a << "\n";
  }
  if (a->bypass_coverage) {
    if (k_cov_debug) {
      std::cout << "Bypass coverage; skipping\n";
    }
    sort_and_output(std::move(a));
    return;
  }
  if (!a->read_coverage) {
    throw(io_exception("pair_cov requires assemblies to be procesesd through read_cov"));
    return;
  }

  track_left_offset(a->left_offset);
  advance_ref_to(a->left_offset);

  if (a->left_offset == a->right_offset) {
    m_cur_inserts.emplace_back(std::move(a));
  } else {
    m_cur_non_inserts.emplace_back(std::move(a));
  }
}

void pair_cov::save_pending_reads(pair_cov_pg* pg) {
  CHECK(pg->pending_pair_table.entries.empty());
  for (const auto& pending : pg->pending_saved_reads) {
    const auto& mm = pending.first;

    pair_table_entry& pte = pg->pending_pair_table.entries[mm];

    for (const auto& pending_pte_entry : pending.second) {
      aoffset_t read_start = pending_pte_entry.first;
      read_id_mask_t mask = pending_pte_entry.second;
      read_id_masks masks(mask, mask);

      result_offset roff;
      roff.r = pg->r.get();
      roff.read_start = read_start;

      pte.results[masks].push_back(roff);
      pte.tot_read_starts +=
          adjust_interval_set(closed_interval(read_start), pg->r->right_anchored_left_offset);
      pte.pending_mask |= mask;
    }
  }

  pg->pending_saved_reads.clear();
}

void pair_cov::save_read(pair_cov_pg* pg, uint32_t read_id, const readmap::read& rd,
                         const multi_mid& mm, aoffset_t read_start) {
  read_id_mask_t read_id_mask = read_id_mask_t(1) << (read_id & (sizeof(read_id_mask_t) * 8 - 1));
  CHECK(read_id_mask);

  pg->pending_saved_reads[mm][read_start] |= read_id_mask;
}

void pair_cov::save_reads(pair_cov_pg* pg) {
  for (const auto& reads_cov_entry : pg->r->a->read_coverage->reads()) {
    aoffset_t read_offset = reads_cov_entry.offset;
    aoffset_t read_len = reads_cov_entry.read_len;
    for (uint32_t read_id : reads_cov_entry.read_ids) {
      readmap::read cov_rd = m_opts.readmap->get_read_by_id(read_id);
      if (cov_rd.is_original_orientation() != m_opts.forward_pairs_face_inward) {
        continue;
      }

      if (!cov_rd.has_mate()) {
        continue;
      }

      // We should see the mate later.
      if (k_cov_debug) {
        std::cout << "Saving read id " << read_id << ": " << cov_rd.get_seqset_entry().sequence()
                  << "\n";
      }

      uint64_t mid_id = cov_rd.get_mid_id();
      CHECK_LT(mid_id, uint64_t(std::numeric_limits<uint32_t>::max()));
      multi_mid mm;
      mm.multi_id = mid_id;
      mm.size = read_len;
      mm.read_id_chunk = read_id / (sizeof(read_id_mask_t) * 8);

      save_read(pg, read_id, cov_rd, mm, read_offset);
    }
  }

  save_pending_reads(pg);
}

void pair_cov::match_mates(pair_cov_pg* pg) {
  // Where this assembly starts in the left-anchored coordinate system.
  aoffset_t seq_start_from_left = pg->r->a->left_offset;
  // Where this assembly starts in the right-anchored coordinate system.
  aoffset_t seq_start_from_right = pg->r->a->right_offset - aoffset_t(pg->r->a->seq.size());

  for (const auto& reads_cov_entry : pg->r->a->read_coverage->reads()) {
    aoffset_t read_offset = reads_cov_entry.offset;
    aoffset_t read_len = reads_cov_entry.read_len;
    for (uint32_t read_id : reads_cov_entry.read_ids) {
      readmap::read cov_rd = m_opts.readmap->get_read_by_id(read_id);
      if (cov_rd.is_original_orientation() == m_opts.forward_pairs_face_inward) {
        continue;
      }
      if (!cov_rd.has_mate()) {
        continue;
      }
      if (k_cov_debug) {
        std::cout << "read id " << read_id
                  << ", looking up read: " << cov_rd.get_seqset_entry().sequence() << "\n";
      }
      // We should have seen the mate for this read previously.
      readmap::read mate = cov_rd.get_mate_rc();
      uint32_t mate_id = mate.get_read_id();
      int mate_size = mate.size();
      uint64_t mid_id = mate.get_mid_id();
      CHECK_LT(mid_id, uint64_t(std::numeric_limits<uint32_t>::max()));
      multi_mid mm;
      mm.multi_id = mid_id;
      mm.size = mate_size;
      mm.read_id_chunk = mate_id / (sizeof(read_id_mask_t) * 8);

      if (k_cov_debug) {
        std::cout << "Mate is: " << mate.get_seqset_entry().sequence() << "\n";
      }

      for (bool pending : {false, true}) {
        pair_table* pt = pending ? &pg->pending_pair_table : &m_main_pair_table;
        aoffset_t seq_start = pending ? seq_start_from_right : seq_start_from_left;
        aoffset_t end_of_read_offset = seq_start + read_offset + read_len;

        auto it = pt->entries.find(mm);
        if (it == pt->entries.end()) {
          if (k_cov_debug) {
            std::cout << "(not present in " << (pending ? "pending" : "main") << ")\n";
          }
          continue;
        }
        auto& pte = it->second;
        if (k_cov_debug) {
          std::cout << "Mate present in" << (pending ? " pending" : " main")
                    << " pair table: " << &pte << "\n";
        }

        if (find_mate(pg, mate_id, mate_size, pte, end_of_read_offset)) {
          // Pair match!
          if (k_cov_debug) {
            std::cout << "Mate found!  Adding read id " << read_id << "\n";
          }
          pg->r->pair_supported_reads.insert(read_offset, read_id, read_len);
        } else {
          if (k_cov_debug) {
            std::cout << "No mate found searching for " << mate_id << "\n";
          }
        }
      }
    }
  }
}

void pair_cov::add_assembly(pair_cov_pg* pg, assembly_ptr a) {
  ++m_cur_asm_count;

  CHECK(pg);
  CHECK(a);

  // There are two coordinate systems.  One where a->left_anchor is
  // anchored to reference, and read offsets are calculated from the
  // beginning, and one where a->right_anchor is anchored to reference
  // and read offsets are calculated from the end.
  //
  // Where svlen is 0, these two coordinate systems are the same.
  //
  // The left anchored coordinate system is used when looking up
  // previously seen mates from assemblies before this.
  //
  // The right anchored coordinate system is used when storing pairs
  // for future use, and for looking up previously seen mates from
  // this assembly.

  // Where this assembly starts in the right-anchored coordinate system.
  aoffset_t seq_start_from_right = a->right_offset - aoffset_t(a->seq.size());

  if (k_cov_debug) {
    std::cout << "\nAdding assembly " << *a << " to table: " << dump(pg);
  }

  CHECK_EQ(m_cur_ref_offset, a->left_offset);

  aoffset_t svlen = aoffset_t(a->seq.size()) - (a->right_offset - a->left_offset);
  if (svlen) {
    aoffset_t care_dist = max_care_distance();
    for (const auto& pend : m_pending_results) {
      aoffset_t dist = (m_cur_ref_offset + aoffset_t(a->seq.size())) -
                       interval_set_upper_bound(pend->cur_sv_adjust);
      CHECK_GE(dist, 0);
      if (dist < care_dist) {
        interval_t adjusted = adjust_interval_set(pend->cur_sv_adjust, closed_interval_set(-svlen));
        auto& adj = m_future_adjusts[a->right_offset][pend.get()];
        if (adj.first) {
          CHECK_EQ(adj.first.get(), pend.get());
          adj.second += adjusted;
        } else {
          adj.first = pend.clone();
          adj.second = adjusted;
        }
      }
    }
  }

  CHECK_EQ(a->left_offset, m_cur_ref_offset);

  CHECK(!pg->r);
  pg->r = result_ptr::make_shared();
  pg->r->left_offset = a->left_offset;
  pg->r->right_offset = a->right_offset;
  pg->r->right_anchored_left_offset = seq_start_from_right;
  pg->r->a = std::move(a);

  save_reads(pg);

  match_mates(pg);

  if (k_cov_debug) {
    std::cout << "Done adding assembly to pair table: " << dump(pg);
  }
}

bool pair_cov::find_mate(pair_cov_pg* pg, uint32_t mate_id, int mate_len, pair_table_entry& pte,
                         aoffset_t end_of_read_offset) {
  interval_set_t pair_valid =
      closed_interval(end_of_read_offset - aoffset_t(m_opts.max_pair_distance),
                      end_of_read_offset - aoffset_t(m_opts.min_pair_distance));

  if (k_cov_debug) {
    std::cout << "find_mate: pair valid interval is " << pair_valid << ", end of read offset is "
              << end_of_read_offset << ", pending mask = " << pte.pending_mask << "\n";
  }

  read_id_mask_t mate_id_mask = read_id_mask_t(1) << (mate_id & (sizeof(read_id_mask_t) * 8 - 1));
  CHECK(mate_id_mask);

  bool found_any = false;
  if (overlaps(pte.tot_read_starts, pair_valid)) {
    found_any = true;
  }

  if (found_any && !(pte.pending_mask & mate_id_mask)) {
    if (k_cov_debug) {
      std::cout << "PTE Pending mask does not contain " << mate_id_mask
                << "; not doing full search\n";
    }
    return found_any;
  }

  bool found_all = true;
  pte_results_t new_results;
  for (auto pm_it = pte.results.begin(); pm_it != pte.results.end();) {
    const auto& masks = pm_it->first;
    auto& roffs = pm_it->second;

    if (masks.total < mate_id_mask) {
      break;
    }

    if (!(masks.total & mate_id_mask)) {
      ++pm_it;
      continue;
    }
    if (found_any && !(masks.pending & mate_id_mask)) {
      ++pm_it;
      continue;
    }

    read_id_mask_t new_pending_mask = masks.pending & ~mate_id_mask;

    size_t idx = 0;
    while (idx < roffs.size()) {
      auto& roff = roffs[idx];

      interval_set_t adjusted_read_start = adjust_interval_set(
          roff.r->cur_sv_adjust, roff.read_start + roff.r->right_anchored_left_offset);
      if (overlaps(pair_valid, adjusted_read_start)) {
        found_any = true;
        pte.tot_read_starts += adjusted_read_start;
        if (masks.pending & mate_id_mask) {
          roff.r->pair_supported_reads.insert(roff.read_start, mate_id, mate_len);
          read_id_masks new_masks(masks.total, new_pending_mask);
          auto existing_it = pte.results.find(new_masks);
          if (existing_it == pte.results.end()) {
            new_results[new_masks].emplace_back(std::move(roff));
          } else {
            existing_it->second.emplace_back(std::move(roff));
          }
          if (idx + 1 != roffs.size()) {
            roffs[idx] = std::move(roffs.back());
          }
          roffs.pop_back();
        } else {
          break;  // Found one, and nothing pending; stop searching.
        }
      } else {
        if (masks.pending & mate_id_mask) {
          found_all = false;
        }
        ++idx;
      }
    }

    if (roffs.empty()) {
      pm_it = pte.results.erase(pm_it);
    } else {
      ++pm_it;
    }
  }

  pte.results.merge(std::move(new_results));

  if (found_all) {
    pte.pending_mask &= ~mate_id_mask;
  }

  return found_any;
}

void pair_cov::reap_result(std::unique_ptr<result> r) {
  if (!r) {
    return;
  }
  if (k_cov_debug) {
    std::cout << "Reaping assembly: " << *r->a << "\n";
  }

  r->a->pair_read_coverage.emplace(r->pair_supported_reads.build_and_clear(r->a->seq.size()));
  untrack_left_offset(r->a->left_offset);
  sort_and_output(std::move(r->a));
}

void pair_cov::join(std::vector<std::unique_ptr<pair_cov_pg>> inputs) {
  for (auto& input : inputs) {
    join(std::move(input));
  }
}

void pair_cov::update_adjusts() {
  while (!m_future_adjusts.empty()) {
    auto it = m_future_adjusts.begin();
    CHECK_GE(it->first, m_cur_ref_offset);
    if (it->first > m_cur_ref_offset) {
      break;
    }

    auto& adjust = it->second;
    for (auto& adj_item : adjust) {
      CHECK_EQ(adj_item.first, adj_item.second.first.get());
      result_ptr r = std::move(adj_item.second.first);
      interval_set_t intv = adj_item.second.second;
      interval_set_t orig_sv_adjust = r->cur_sv_adjust;
      r->cur_sv_adjust += intv;
      // We shouldn't adjust the valid placement range of this old
      // assembly to be past our current position.
      aoffset_t adjust_limit = m_cur_ref_offset - r->right_offset;

      aoffset_t adjust_upper_closed = r->cur_sv_adjust.upper() - 1;
      CHECK_LE(adjust_upper_closed, adjust_limit)
          << " updating pending result " << *r->a << ": orig = " << orig_sv_adjust
          << ", intv = " << intv << ", after adjust = " << r->cur_sv_adjust << "\n";

      if (!m_pending_results.count(r)) {
        m_pending_results.insert(std::move(r));
      } else {
        CHECK(!r.release());
      }
    }

    m_future_adjusts.erase(it);
  }
}

void pair_cov::join(std::unique_ptr<pair_cov_pg> input) {
  propagate_and_fill(input->pending_pair_table, &m_main_pair_table);
  m_pending_results.insert(std::move(input->r));
}

pair_cov::simple_interval_set pair_cov::adjust_interval_set(const simple_interval_set& orig,
                                                            aoffset_t adjust) {
  return simple_interval_set(orig.lower() + adjust, orig.upper() + adjust);
}

pair_cov::simple_interval_set pair_cov::adjust_interval_set(const simple_interval_set& orig,
                                                            const simple_interval_set& adjust) {
  aoffset_t adjust_upper_max = adjust.upper() - 1;
  return simple_interval_set(orig.lower() + adjust.lower(), orig.upper() + adjust_upper_max);
}

pair_cov::simple_interval_set pair_cov::unadjust_interval_set(const simple_interval_set& orig,
                                                              const simple_interval_set& adjust) {
  aoffset_t adjust_upper_max = adjust.upper() - 1;
  return simple_interval_set(orig.lower() - adjust_upper_max, orig.upper() - adjust.lower());
}

aoffset_t pair_cov::max_care_distance() const {
  aoffset_t read_dist = aoffset_t(m_opts.seqset->max_read_len() * 2);
  aoffset_t pair_dist = aoffset_t(m_opts.max_pair_distance);
  return read_dist + pair_dist;
}
absl::flat_hash_set<pair_cov::result*> pair_cov::find_expired_results() const {
  aoffset_t pair_cutoff = m_cur_ref_offset - max_care_distance();

  if (k_cov_debug) {
    std::cout << "Finding expired results, ref=" << m_cur_ref_offset
              << " pair cutoff=" << pair_cutoff << "\n";
  }

  absl::flat_hash_set<result*> expired;
  for (const auto& r : m_pending_results) {
    aoffset_t adjusted_offset_from_cur = interval_set_upper_bound(
        adjust_interval_set(closed_interval_set(r->right_offset), r->cur_sv_adjust));

    if (k_cov_debug) {
      std::cout << "@" << m_cur_ref_offset << ": right offset: " << r->right_offset
                << " adjusted by " << r->cur_sv_adjust << ": " << adjusted_offset_from_cur;
      if (adjusted_offset_from_cur < pair_cutoff) {
        std::cout << " EXPIRED\n";
      } else {
        std::cout << " OK\n";
      }
    }
    if (adjusted_offset_from_cur < pair_cutoff) {
      expired.insert(r.get());
    }
  }
  return expired;
}

namespace {

template <typename T>
class track_histo {
 public:
  track_histo() = default;

  void increment(const T& val) { ++m_counts[val]; }

  void display_brief_pct(std::ostream& os, size_t display_count) {
    size_t tot = 0;
    for (const auto& val : m_counts) {
      tot += val.second;
    }

    size_t tot_so_far = 0;
    size_t i = 0;
    for (const auto& val : m_counts) {
      if (!display_count) {
        os << "...";
        return;
      }
      float pct_here = val.second * 100. / tot;
      tot_so_far += val.second;
      float pct_so_far = tot_so_far * 100. / tot;
      if (i == 0 || !(i & (i - 1))) {
        os << std::fixed << std::setprecision(2) << " " << val.first << "(" << val.second << "+"
           << pct_here << "%=" << pct_so_far << "%)";
      }
      --display_count;
      ++i;
    }
  }

 private:
  std::map<T, size_t> m_counts;
};

}  // namespace

void pair_cov::flush_old() {
  absl::flat_hash_set<result*> expired_results = find_expired_results();
  size_t expired_entries = 0;
  size_t expired_roff = 0;
  size_t remaining_roff = 0;
  size_t pending_mates = 0;
  size_t expired_pending_mates = 0;

  track_histo<size_t> pte_results_histo;
  track_histo<size_t> roff_per_mask_histo;

  multi_mid most_roffs_pte_mm{};
  pair_table_entry* most_roffs_pte = nullptr;
  size_t most_roffs_pte_pms = 0;

  multi_mid most_results_pte_mm{};
  pair_table_entry* most_results_pte = nullptr;
  size_t most_results_pte_pms = 0;

  auto& pt = m_main_pair_table.entries;
  if (expired_results.empty() && !k_extended_stats) {
    return;
  }

  if (k_cov_debug || k_extended_stats) {
    std::cout << "flush_old, expiring " << expired_results.size() << " from pair table of size "
              << pt.size() << ": ";
    std::cout.flush();
  }
  auto it = pt.begin();
  for (auto next = it; it != pt.end(); it = next) {
    ++next;
    const auto& mm = it->first;
    auto& pte = it->second;

    CHECK(!(pte.results.empty()));
    size_t tot_roffs = 0;

    auto results_it = pte.results.begin();
    bool do_shrink_results = false;
    while (results_it != pte.results.end()) {
      auto& roffs = results_it->second;
      if (k_extended_stats) {
        roff_per_mask_histo.increment(roffs.size());
      }
      bool do_shrink_roffs = false;
      size_t idx = 0;
      tot_roffs += roffs.size();
      while (idx != roffs.size()) {
        auto& roff = roffs[idx];
        if (expired_results.count(roff.r)) {
          do_shrink_roffs = true;
          ++expired_roff;
          if (idx + 1 != roffs.size()) {
            roff = std::move(roffs.back());
          }
          roffs.pop_back();
        } else {
          ++idx;
          ++remaining_roff;
        }
      }
      if (roffs.empty()) {
        results_it = pte.results.erase(results_it);
        do_shrink_results = true;
      } else {
        ++results_it;
        if (do_shrink_roffs) {
          roffs.shrink_to_fit();
        }
      }
    }

    if (pte.results.empty()) {
      pt.erase(it);
      ++expired_entries;
      continue;
    }

    if (do_shrink_results) {
      pte.results.shrink_to_fit();
    }

    if (k_extended_stats) {
      if (!most_roffs_pte || tot_roffs > most_roffs_pte_pms) {
        most_roffs_pte = &pte;
        most_roffs_pte_mm = mm;
        most_roffs_pte_pms = tot_roffs;
      }

      if (!most_results_pte || pte.results.size() > most_results_pte->results.size()) {
        most_results_pte = &pte;
        most_results_pte_mm = mm;
        most_results_pte_pms = tot_roffs;
      }

      pte_results_histo.increment(pte.results.size());
    }
  }

  if (k_cov_debug || k_extended_stats) {
    std::cout << " expired " << expired_results.size() << " results, " << expired_entries
              << " pte, " << expired_roff << "/" << remaining_roff << " roffs, "
              << expired_pending_mates << "/" << pending_mates
              << " pending mates, new size=" << pt.size() << "\n";
    std::cout << "roffs per mask histo: ";
    roff_per_mask_histo.display_brief_pct(std::cout, 64);
    std::cout << "\nresults histo: ";
    pte_results_histo.display_brief_pct(std::cout, 64);
    std::cout << "\n";

    if (most_roffs_pte) {
      std::cout << "Most pms pte has " << most_roffs_pte_pms << " pms and "
                << most_roffs_pte->results.size()
                << " results: " << dump_multi_mid(most_roffs_pte_mm) << "\n";
    }
    if (most_results_pte) {
      std::cout << "Most results pte has " << most_results_pte_pms << " pms and "
                << most_results_pte->results.size()
                << " results: " << dump_multi_mid(most_results_pte_mm) << "\n";
    }
  }

  for (const auto& expired : expired_results) {
    auto it = m_pending_results.find<const result*>(expired);
    CHECK(it != m_pending_results.end());
    result_ptr r = std::move(m_pending_results.extract(it).value());
    reap_result(r.release());
  }

  if (k_extended_stats) {
    result* oldest = nullptr;
    for (const auto& r : m_pending_results) {
      if (!oldest || r->a->left_offset < oldest->a->left_offset) {
        oldest = r.get();
      }
    }

    if (oldest) {
      std::cout << "pair_cov's oldest assembly is at " << oldest->a->left_offset << ", ["
                << m_cur_ref_offset - oldest->a->left_offset << ","
                << m_cur_ref_offset - oldest->a->right_offset << ") behind cur(" << m_cur_ref_offset
                << "), valid adjust range=" << oldest->cur_sv_adjust << " : " << *oldest->a << "\n";
    }

    std::cout << sorted_output_stats(m_cur_ref_offset) << "\n";
  }
}

std::string pair_cov::dump_multi_mid(const multi_mid& mm) const {
  std::stringstream result;
  uint64_t seqset_id = m_opts.readmap->mid_to_entry(mm.multi_id);
  seqset_range r = m_opts.seqset->ctx_entry(seqset_id).truncate(mm.size);

  auto read_range = m_opts.readmap->entry_to_index_range(r.begin(), r.end());
  size_t read_count = read_range.second - read_range.first;
  result << "Mid(id=" << mm.multi_id << ",size=" << mm.size << ",nr=" << read_count
         << ",seq=" << r.sequence() << ")";
  return result.str();
}

void pair_cov::propagate_and_fill(pair_table& old_table, pair_table* new_table) {
  if (k_cov_debug) {
    std::cout << "Propagating and filling from table with " << old_table.entries.size()
              << " entries to table with " << new_table->entries.size() << " entries\n";
  }

  for (auto& old_pte_entry : old_table.entries) {
    const multi_mid& mm = old_pte_entry.first;
    pair_table_entry& old_pte = old_pte_entry.second;

    if (k_cov_debug > 1) {
      std::cout << "Read mid=" << dump_multi_mid(mm) << " old pte=" << &old_pte << "\n";
    }

    pair_table_entry& new_pte = new_table->entries[mm];

    if (k_cov_debug > 1) {
      std::cout << "New pte: " << &new_pte << "\n";
    }

    static size_t count = 0;
    if (k_cov_debug && !(count++ & 0xFFFFF)) {
      std::cout << " results size: " << old_pte.results.size() << "\n";
    }

    new_pte.results.reserve(new_pte.results.size() + old_pte.results.size());

    for (auto& old_pending : old_pte.results) {
      auto& old_roffs = old_pending.second;
      auto& new_roffs = new_pte.results[old_pending.first];
      new_roffs.insert(new_roffs.end(), old_roffs.begin(), old_roffs.end());
    }

    new_pte.tot_read_starts += old_pte.tot_read_starts;
    new_pte.pending_mask |= old_pte.pending_mask;
  }
  if (k_cov_debug) {
    std::cout << "Done propagating and filling\n";
  }
}

pair_cov::interval_t pair_cov::closed_interval(aoffset_t start, aoffset_t limit) {
  return interval_t(start, limit + 1);
}

pair_cov::interval_t pair_cov::closed_interval(aoffset_t single_offset) {
  return closed_interval(single_offset, single_offset);
}

aoffset_t pair_cov::interval_set_upper_bound(const interval_set_t& s) {
  CHECK(!is_empty(s));
  return s.upper();
}

pair_cov::interval_set_t pair_cov::closed_interval_set(aoffset_t single_offset) {
  interval_set_t result;
  result += closed_interval(single_offset);
  return result;
}

}  // namespace variants

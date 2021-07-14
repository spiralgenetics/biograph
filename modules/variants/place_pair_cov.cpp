#include "modules/variants/place_pair_cov.h"
#include "modules/variants/apply_edges.h"

#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"

#include <dlfcn.h>

static constexpr bool k_dbg = false;
static constexpr bool k_dbg_some = false;
static constexpr bool k_stats = false;
static constexpr bool k_dist_stats = false;

namespace variants {

absl::btree_set<uint32_t> place_pair_cov::g_debug_read_ids = {
    /*38852481, 1408171319, 446049297,
                                                                 220893869,
                                                                 //
                                                                 125651015*/
};

absl::btree_set<size_t> place_pair_cov::g_debug_assembly_ids = {};

void place_pair_cov::debug_assembly(const assembly& a) {
  std::cerr << "Debugging reads in assembly " << a << ", read ids:";
  g_debug_assembly_ids.insert(a.assembly_id);
  for (uint32_t read_id : a.read_coverage->all_read_ids()) {
    auto rd = m_opts.readmap->get_read_by_id(read_id);
    if (rd.has_mate()) {
      g_debug_read_ids.insert(read_id);
      std::cerr << " " << read_id;
    } else {
      std::cerr << " (mateless:" << read_id << ")";
    }
  }
  std::cerr << "\n";
}

void place_pair_cov::maybe_debug_assembly(const assembly& a) {
  if (!k_dbg_some) {
    return;
  }
  aoffset_t left = a.left_offset;
  aoffset_t right = a.right_offset;

  //  if (left == 669651 && right == 669651) {
  //    // chr1, Insert length 231 should be the proper one.
  //    debug_assembly(a);
  if (left >= 25528100 && right <= 25528125 &&
      (a.seq ==
           dna_sequence("GGACAGTGATGGCGGAGGTGGTGGAGTCGGTGGAAGGACAGTGATGGCGGAGGTTGTGGAGTTGGTGGAGGGAC"
                        "AGTGATGGTGGAGGTGGTGGAGTCGGTGGAGGGACAGTGATGGTGGAGGTGGTGGAGTCGGTGGATGGACAGTG"
                        "CGGGTGGAGGTAGTGGAGTCGGTGGAGGGACAGTGATGGTGGAGGTGGTGGAGTCGGTGGAGGGACAGTGATGG"
                        "TGGAGGTGGTGGAGTCGGTGGAGGGACAGTGATGGTGGAGGTGGTGGAGTCGGTGGAG") ||
       a.seq == dna_sequence("GGTGGAGGTGGTGGAGTCGGTGGATGGACAGTGCGGGTGGAGGTAGTGGAGTTGGTGGAGGGACAGTGA"
                             "TGGTGGAGGTGGTGGAGTCGGTGGAGGGACAGTGATGGTGGAGGTGGTGGAGTCGGTGGAGGGACAGTG"
                             "ATGGTGGAGGTGGTGGAGTCGGTGGAGGGACAGTGCGGGTGGAGGTGGTGGAGTCGGTGGAGGGACAGT"
                             "GATGGTGGAGGTGGTGGAGTCGGTGGAGGGACAGTGCG"))) {
    // chr15	25528119	25528120	REF	INS	[200,300)
    debug_assembly(a);
  } else {
    if (!a.matches_reference) {
      std::cerr << "Not debugging " << a << "\n";
    }
  }
}

place_pair_cov::place_pair_cov(const assemble_options& opts, const place_pair_options& popts,
                               pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_opts(opts), m_popts(popts) {
  aoffset_t min_read_len = m_opts.readmap->min_read_len();
  aoffset_t max_read_len = m_opts.readmap->max_read_len();
  m_max_dist = std::max<aoffset_t>(aoffset_t(m_opts.max_pair_distance) - 2 * min_read_len,
                                   // Overlapping:
                                   max_read_len);

  m_max_ideal_dist = std::max<aoffset_t>(aoffset_t(m_popts.ideal_pair_distance) - 2 * min_read_len,
                                         // Overlapping:
                                         max_read_len);
}

place_pair_cov::~place_pair_cov() { CHECK(m_block.empty()); }

void place_pair_cov::flush() {
  if (m_block.empty()) {
    return;
  }
  place();
  if (m_inspect_for_testing) {
    m_inspect_for_testing(this);
  }
  for (auto& a : m_block) {
    sort_and_output(std::move(a));
  }
  m_block.clear();
  flush_sorted();
}

void place_pair_cov::place() {
  aoffset_t last_pos = 0;
  size_t last_bytes = 0;
  const char k_cur_alloc_prop_name[] = "generic.current_allocated_bytes";
  typedef int (*GetNumericPropertyFunc)(const char*, size_t*);
  GetNumericPropertyFunc GetMallocProperty = nullptr;

  if (k_stats) {
    if (!m_block.empty()) {
      for (const auto& a : m_block) {
        last_pos = std::max<aoffset_t>(last_pos, a->right_offset);
      }
    }

    // Retrieve current allocation total if possible.
    GetMallocProperty =
        ((GetNumericPropertyFunc)dlsym(RTLD_DEFAULT, "MallocExtension_GetNumericProperty"));
    if (GetMallocProperty) {
      CHECK(GetMallocProperty(k_cur_alloc_prop_name, &last_bytes));
    }
  }
  auto last_time = std::chrono::high_resolution_clock::now();
  auto timediff = [&]() -> std::string {
    auto cur_time = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds ms_dist =
        std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - last_time);
    last_time = cur_time;

    size_t cur_bytes = 0;
    if (GetMallocProperty) {
      CHECK(GetMallocProperty(k_cur_alloc_prop_name, &cur_bytes));
    }
    ssize_t bytes_diff = ssize_t(cur_bytes) - ssize_t(last_bytes);
    last_bytes = cur_bytes;

    return absl::StrFormat("%.2fs %+.2fMB", ms_dist.count() / 1000., bytes_diff / (1024 * 1024.));
  };

  if (k_dbg_some) {
    for (const auto& a : m_block) {
      maybe_debug_assembly(*a);
    }
    if (!g_debug_read_ids.empty()) {
      absl::btree_set<uint32_t> ids = g_debug_read_ids;
      for (uint32_t read_id : ids) {
        auto rd = m_opts.readmap->get_read_by_id(read_id);
        CHECK(rd.has_mate()) << read_id;
        g_debug_read_ids.insert(rd.get_rev_comp().get_read_id());
        g_debug_read_ids.insert(rd.get_mate_rc().get_read_id());
        g_debug_read_ids.insert(rd.get_mate().get_read_id());
      }
    }
  }
  dump_state("place start");
  if (k_stats) {
    std::cerr << "Start " << last_pos << "/";
  }
  calc_dists();
  dump_state("calc dists");
  if (k_stats) {
    std::cerr << "CalcDist " << last_pos << " " << timediff() << "/";
  }
  init_edges();
  dump_state("init edges");
  if (k_stats) {
    std::cerr << "InitEdge " << last_pos << " " << timediff() << "/";
  }
  filter_reads();
  dump_state("filter reads");
  if (k_stats) {
    std::cerr << "FiltReads " << last_pos << " " << timediff() << "/";
  }
  save_filtered_reads();
  if (k_stats) {
    std::cerr << "DoneFilt " << last_pos << " " << timediff() << "\n";
  }
}

place_pair_cov::trace_state* place_pair_cov::assembly_to_state(assembly* a) {
  auto it = m_asm_to_state.find(a);
  CHECK(it != m_asm_to_state.end());
  return &*it->second;
}

void place_pair_cov::add_dists(const dists_t& old, aoffset_t distance, dists_t& result) {
  CHECK_GE(distance, 0);
  if (k_dbg) {
    std::cerr << "Adding " << distance << " distance to " << dump_dist_table(old) << " into "
              << dump_dist_table(result) << ", max = " << m_max_dist
              << " max ideal=" << m_max_ideal_dist << "\n";
  }
  for (const auto& old_elem : old) {
    const auto& old_offset = old_elem.first;
    const auto& old_dists = old_elem.second;

    dist_set added = old_dists.add_offset(distance, m_max_dist, m_max_ideal_dist);
    if (!added.empty()) {
      result[old_offset] |= added;
    }
  }
  if (k_dbg) {
    std::cerr << "New dist table is  " << dump_dist_table(result) << "\n";
  }
}

void place_pair_cov::calc_dists() {
  aoffset_t last_offset = std::numeric_limits<aoffset_t>::min();
  apply_edges_to_block(  //
      m_block, [this, &last_offset](aoffset_t offset, const std::vector<assembly_ptr>& left_edges,
                                    const std::vector<assembly_ptr>& inserts,
                                    const std::vector<assembly_ptr>& right_edges) {
        auto& cur_dists = m_dists[offset];

        // Dist to self is 0
        cur_dists[offset].insert(0);

        // Always add ref distance
        CHECK_GT(offset, last_offset);
        if (last_offset != std::numeric_limits<aoffset_t>::min()) {
          auto it = m_dists.find(last_offset);
          CHECK(it != m_dists.end());

          add_dists(it->second, (offset - last_offset), cur_dists);
        }
        last_offset = offset;

        for (const auto& left : left_edges) {
          auto it = m_dists.find(left->left_offset);
          CHECK(it != m_dists.end());

          add_dists(it->second, left->seq.size(), cur_dists);
        }

        if (!inserts.empty()) {
          dists_t after_inserts;
          for (const auto& insert : inserts) {
            add_dists(cur_dists, insert->seq.size(), after_inserts);
          }
          add_dists(after_inserts, 0, cur_dists);
        }

        if (k_dist_stats && !cur_dists.empty()) {
          size_t tot = 0;
          aoffset_t tot_min = std::numeric_limits<aoffset_t>::max();
          aoffset_t tot_max = std::numeric_limits<aoffset_t>::min();

          for (const auto& d : cur_dists) {
            CHECK(!d.second.empty());
            tot += d.second.size();
            tot_min = std::min(tot_min, d.second.min_value());
            tot_max = std::max(tot_max, d.second.max_value());
          }

          std::cerr << "@" << offset << ": Cur dists " << cur_dists.size() << " ("
                    << cur_dists.begin()->first << "," << cur_dists.rbegin()->first
                    << "), tot=" << tot << " (" << tot_min << "," << tot_max << ")\n";
        }
      });
}

void place_pair_cov::save_reads(assembly* a, trace_state* st) {
  auto* rm = m_opts.readmap;
  for (const auto& cov_entry : a->read_coverage->reads()) {
    // Find left mates that end here
    if (cov_entry.offset + cov_entry.read_len <= aoffset_t(a->seq.size())) {
      for (uint32_t read_id : cov_entry.read_ids) {
        if (!rm->has_mate(read_id)) {
          continue;
        }
        if (rm->get_is_forward(read_id) != m_opts.forward_pairs_face_inward) {
          continue;
        }
        anchor left;
        left.st = st;
        left.offset = cov_entry.offset;
        m_reads[read_id].read_states.emplace(left);
      }
    }
    // Find right mates that start here
    if (cov_entry.offset >= 0) {
      for (uint32_t read_id : cov_entry.read_ids) {
        if (!rm->has_mate(read_id)) {
          continue;
        }
        if (rm->get_is_forward(read_id) == m_opts.forward_pairs_face_inward) {
          continue;
        }
        anchor right;
        right.st = st;
        right.offset = cov_entry.offset;
        uint32_t left_mate_id = rm->get_mate_rc(read_id);
        m_reads[left_mate_id].rc_mate_states.emplace(right);
      }
    }
  }
}

void place_pair_cov::init_edges() {
  apply_edges_to_block(  //
      m_block, [this](aoffset_t, const std::vector<assembly_ptr>& left_edges,
                      const std::vector<assembly_ptr>& inserts,
                      const std::vector<assembly_ptr>& right_edges) {
        auto lookup_and_add = [this](const std::vector<assembly_ptr>& edges,
                                     std::vector<trace_state*>& dest) {
          for (const auto& a : edges) {
            dest.push_back(assembly_to_state(a.get()));
          }
        };
        for (auto& a : left_edges) {
          trace_state* state = assembly_to_state(a.get());
          lookup_and_add(inserts, state->right_edges);
          lookup_and_add(right_edges, state->right_edges);
        }
        for (auto& a : inserts) {
          trace_state* state = assembly_to_state(a.get());
          lookup_and_add(left_edges, state->left_edges);
          lookup_and_add(right_edges, state->right_edges);
        }
        for (auto& a : right_edges) {
          trace_state* state = assembly_to_state(a.get());
          lookup_and_add(left_edges, state->left_edges);
          lookup_and_add(inserts, state->left_edges);
        }
      });
}

void place_pair_cov::filter_reads() {
  size_t tot_reads = m_reads.size();
  size_t reads_so_far = 0;
  size_t report_every = std::max<size_t>(tot_reads / 10, 1);
  for (auto& elem : m_reads) {
    uint32_t read_id = elem.first;
    read_info ri = std::move(elem.second);
    place_and_filter(read_id, ri);
    if (k_stats) {
      ++reads_so_far;
      if ((reads_so_far % report_every) == 0) {
        std::cerr << "Placing read id=" << read_id << ", #" << reads_so_far << "/" << tot_reads
                  << " (" << (reads_so_far * 100 / tot_reads) << "%\n";
      }
    }
  }
  m_reads.clear();
}

void place_pair_cov::gather_anchor(gather_anchors& anchors, const pair_align_metric& metric,
                                   const anchor* left, const anchor* right, bool brief_dbg) {
  bool dbg = k_dbg;

  if (metric.dist_from_ideal > (m_opts.max_pair_distance - m_popts.ideal_pair_distance) ||
      metric.dist_from_ideal < (m_opts.min_pair_distance - m_popts.ideal_pair_distance)) {
    // Enforce maximum and minimum bounds.
    if (dbg) {
      std::cerr << "Metric out of range: " << metric << "\n";
    }
    if (brief_dbg) {
      std::cerr << "R";
    }
    return;
  }
  if (dbg) {
    std::cerr << "Checking metric: " << metric << "\n";
  }
  if (metric_better(anchors.best_metric, metric)) {
    if (brief_dbg) {
      std::cerr << " WORSE";
    }
    if (dbg) {
      std::cerr << "Old best is better: " << anchors.best_metric << "\n";
    }
    CHECK(!anchors.best_pairs.empty());
    return;
  }

  if (metric_better(metric, anchors.best_metric)) {
    if (brief_dbg) {
      std::cerr << " NEW";
    }
    if (dbg) {
      std::cerr << "New best is better: " << metric << "\n";
    }
    anchors.best_pairs.clear();
    anchors.best_metric = metric;
  }
  if (brief_dbg) {
    std::cerr << " BEST";
  }
  anchors.best_pairs.emplace_back(left, right);
  if (brief_dbg) {
    std::cerr << "\n";
  }
}

void place_pair_cov::place_and_filter(uint32_t read_id, const read_info& ri) {
  bool dbg = k_dbg;
  bool brief_dbg = k_dbg;

  //  if (read_id == 44297) {
  //    dbg = true;
  //  }

  if (g_debug_read_ids.count(read_id)) {
    brief_dbg = true;
  }

  if (dbg) {
    std::cerr << "Placing and filtering read " << read_id << "\n";
  }

  m_lefts.clear();
  m_rights.clear();
  auto* rm = m_opts.readmap;
  uint32_t rc_mate_id = rm->get_mate_rc(read_id);
  aoffset_t read_len = rm->get_readlength(read_id);
  aoffset_t mate_len = rm->get_readlength(rc_mate_id);

  const auto& lefts = ri.read_states;
  const auto& rights = ri.rc_mate_states;

  if (lefts.empty() || rights.empty()) {
    if (dbg) {
      std::cerr << "Missing mates for read " << read_id << " with " << lefts.size() << " lefts and "
                << rights.size() << " rights\n";
    }
    return;
  }

  if (brief_dbg) {
    std::cerr << "Read " << read_id << " Total: " << lefts.size() << " left and " << rights.size()
              << " rights\n";
  }

  if (dbg) {
    for (const auto& right : rights) {
      std::cerr << "Right at " << right.offset << ": " << *right.st << "\n";
    }
    for (const auto& left : lefts) {
      std::cerr << "Left at " << left.offset << ": " << *left.st << "\n";
    }
  }
  if (dbg) {
    std::cerr << "Left read: " << rm->get_read_by_id(read_id).get_seqset_entry().sequence() << "\n";
    std::cerr << "Right read: " << rm->get_read_by_id(rc_mate_id).get_seqset_entry().sequence()
              << "\n";
  }

  gather_anchors anchors;
  for (const auto& right : rights) {
    dists_table_t::iterator right_dist_it = m_dists.end();
    boost::optional<absl::btree_set<aoffset_t /* overlap bases between mates */>> right_overlaps;

    auto* right_a = right.st->a;
    for (const auto& left : lefts) {
      auto* left_a = left.st->a;

      if (dbg) {
        std::cerr << "Considering anchoring with:\nLeft:  " << *left.st << "\nRight: " << *right.st
                  << "\n";
      }

      if (left.st == right.st) {
        if (dbg) {
          std::cerr << "Synthesizing same-assembly ref distance: " << left.st->a->seq.size()
                    << "\n";
        }
        // Distance from left read end to the end of the assembly
        aoffset_t left_dist = aoffset_t(left_a->seq.size()) - (left.offset + read_len);
        // Distance from beginning of the assembly to the right read start.
        aoffset_t right_dist = right.offset;

        int ref_dist = -left_a->seq.size();
        pair_align_metric metric;
        aoffset_t relative_ideal =
            m_popts.ideal_pair_distance - (read_len + left_dist + right_dist + mate_len);
        metric.dist_from_ideal = ref_dist - relative_ideal;
        gather_anchor(anchors, metric, &left, &right, brief_dbg);
        continue;
      } else if (left_a->right_offset <= right_a->left_offset) {
        // Non-overlapping case; read ends in an assembly that's before mate.

        if (right_dist_it == m_dists.end()) {
          right_dist_it = m_dists.find(right_a->left_offset);
          CHECK(right_dist_it != m_dists.end());
        }

        const auto& right_dists = right_dist_it->second;
        auto dist_it = right_dists.find(left_a->right_offset);
        if (dist_it == right_dists.end()) {
          if (dbg) {
            std::cerr << "No distances available between " << right_a->left_offset << " and "
                      << left_a->right_offset << ", left_st = " << *left.st
                      << ", right_st = " << *right.st << "\n";
          }
          CHECK_NE(left_a->right_offset, right_a->left_offset);
          if (brief_dbg) {
            if (left_a->right_offset < right_a->left_offset) {
              std::cerr << " NODIST, too far ";
            } else {
              std::cerr << " NODIST, wrong side ";
            }
            std::cerr << right_a->left_offset - left_a->right_offset << "\n";
          }
          continue;
        }
        const auto& ref_dists = dist_it->second;
        if (ref_dists.empty()) {
          if (brief_dbg) {
            std::cerr << " NODIST, empty\n";
          }
          continue;
        }

        // Distance from left read end to the end of the assembly
        aoffset_t left_dist = aoffset_t(left_a->seq.size()) - (left.offset + read_len);
        // Distance from beginning of the assembly to the right read start.
        aoffset_t right_dist = right.offset;

        CHECK_GE(left_dist, 0);
        CHECK_LT(left_dist, aoffset_t(left_a->seq.size()));

        CHECK_GE(right_dist, 0);
        CHECK_LT(right_dist, aoffset_t(right_a->seq.size()));

        aoffset_t relative_ideal =
            m_popts.ideal_pair_distance - (read_len + left_dist + right_dist + mate_len);

        int ref_dist = ref_dists.closest_distance_to(relative_ideal);
        pair_align_metric metric;
        metric.dist_from_ideal = ref_dist - relative_ideal;

        if (dbg) {
          std::cerr << "Non-overlapping case, left dist = " << left_dist
                    << ", right_dist = " << right_dist << ", ref_dist = " << ref_dist
                    << " relative ideal= " << relative_ideal << "\n";
        }
        if (brief_dbg) {
          std::cerr << ", " << read_len + left_dist + ref_dist + right_dist + mate_len << ":";
        }
        gather_anchor(anchors, metric, &left, &right, brief_dbg);
      } else if (right_a->right_offset <= left_a->left_offset) {
        // Overlapping case, where the read ends after the mate
        // begins.  The region between the two assemblies is present
        // in both reads.

        // Distance from right read start to end of assembly
        aoffset_t right_dist = aoffset_t(right_a->seq.size()) - right.offset;
        CHECK_GT(right_dist, 0);

        if (!right_overlaps) {
          right_overlaps.emplace();
          aoffset_t max_allowed_overlap =
              (read_len + mate_len) - aoffset_t(m_opts.min_pair_distance);
          if (dbg) {
            std::cerr << "Maximum allowed overlap: " << max_allowed_overlap << "\n";
          }
          if (max_allowed_overlap <= 0) {
            continue;
          }
          for (const auto& cov_entry : right_a->read_coverage->reads()) {
            // Left read must end in a different assembly to the right; otherwise
            // we would catch this case with the same-assembly case.
            if (cov_entry.offset + read_len < aoffset_t(right_a->seq.size())) {
              continue;
            }
            aoffset_t overlap = cov_entry.offset + read_len - right.offset;
            if (overlap <= right_dist || overlap > max_allowed_overlap) {
              continue;
            }

            // Check if the left read ends after the right read starts.
            if (cov_entry.read_ids.contains(read_id)) {
              right_overlaps->insert(overlap);
            }
          }
          if (dbg) {
            std::cerr << "Generated right overlap table: " << StrJoin(*right_overlaps, ",") << "\n";
          }
        }

        if (right_overlaps->empty()) {
          if (dbg) {
            std::cerr << "No right overlaps available.\n";
          }
          continue;
        }

        auto left_ol_dist_it = m_dists.find(left_a->left_offset);
        CHECK(left_ol_dist_it != m_dists.end());

        const auto& left_dists = left_ol_dist_it->second;

        auto dist_it = left_dists.find(right_a->right_offset);
        if (dist_it == left_dists.end()) {
          if (dbg) {
            std::cerr << "No overlapping distances available between " << left_a->left_offset
                      << " and " << right_a->right_offset << ", left_st = " << *left.st
                      << ", right_st = " << *right.st << "\n";
          }
          continue;
        }
        const auto& ref_dists = dist_it->second;

        if (ref_dists.empty()) {
          if (brief_dbg) {
            std::cerr << " NODIST, empty\n";
          }
          continue;
        }
        // Distance from beginning of assembly to left read end
        aoffset_t left_dist = left.offset + read_len;
        CHECK_GT(left_dist, 0);

        for (aoffset_t overlap : *right_overlaps) {
          aoffset_t ref_dist = overlap - (right_dist + left_dist);
          if (dbg) {
            std::cerr << "Trying overlap " << overlap << ", ref dist " << ref_dist
                      << ", right dist " << right_dist << ", left dist " << left_dist << "\n";
          }
          if (!ref_dists.contains(ref_dist)) {
            if (dbg) {
              std::cerr << "Ref dist  " << ref_dist << " not contained in ref dists " << ref_dists
                        << "\n";
            }
            continue;
          }
          if (!left_a->read_coverage->get_read_ids_at(left.offset + read_len - overlap, mate_len)
                   .contains(rc_mate_id)) {
            if (dbg) {
              std::cerr << "Left doesn't support this overlap\n";
            }
            continue;
          }
          pair_align_metric metric;
          aoffset_t pair_dist = read_len + mate_len - overlap;

          metric.dist_from_ideal = pair_dist - m_popts.ideal_pair_distance;

          if (dbg) {
            std::cerr << "Overlapping case, overlap = " << overlap
                      << ", dist from ideal = " << metric.dist_from_ideal << "\n";
          }
          gather_anchor(anchors, metric, &left, &right, brief_dbg);
        }
      } else {
        if (dbg) {
          std::cerr << "Neither overlapping or non-overlapping\n";
        }
      }
    }
  }

  if (anchors.best_pairs.empty()) {
    // No matches found.
    if (dbg || brief_dbg) {
      std::cerr << "No matches found\n";
    }
    return;
  }

  if (anchors.best_pairs.size() > m_popts.max_ambig) {
    if (dbg || brief_dbg) {
      std::cerr << "Skipping ambiguous; best size=" << anchors.best_pairs.size()
                << ", max=" << m_popts.max_ambig << "\n";
    }
    return;
  }

  std::vector<std::pair<const align*, const align*>> best_aligns;

  // Save align pointers so we don't have to copy align structures.
  std::vector<std::vector<align>> align_storage;

  for (const auto& best : anchors.best_pairs) {
    const anchor* left_anchor = best.first;
    const anchor* right_anchor = best.second;

    std::vector<align> left_aligns;
    CHECK(m_so_far.parts.empty());
    bool did_propagate = propagate_read(read_id, left_anchor->st, left_anchor->offset, read_len,
                                        true /* propagate left */, left_aligns, m_so_far);
    CHECK(did_propagate);
    CHECK(m_so_far.parts.empty());

    std::vector<align> right_aligns;
    did_propagate = propagate_read(rc_mate_id, right_anchor->st, right_anchor->offset, mate_len,
                                   false /* propagate right */, right_aligns, m_so_far);
    CHECK(did_propagate);
    CHECK(m_so_far.parts.empty());

    for (const auto& left_align : left_aligns) {
      for (const auto& right_align : right_aligns) {
        best_aligns.emplace_back(&left_align, &right_align);
        if (best_aligns.size() > m_popts.max_ambig) {
          if (dbg || brief_dbg) {
            std::cerr << "Skipping ambiguous after align; best size=" << anchors.best_pairs.size()
                      << " aligns=" << best_aligns.size() << ", max=" << m_popts.max_ambig << "\n";
          }
          return;
        }
      }
    }

    align_storage.emplace_back(std::move(left_aligns));
    align_storage.emplace_back(std::move(right_aligns));
  }

  if (k_dbg) {
    std::cerr << "With read id " << read_id << ", best_metric = " << anchors.best_metric
              << ", bounds = "
              << (aoffset_t(m_opts.max_pair_distance) - m_popts.ideal_pair_distance) << " to "
              << (aoffset_t(m_opts.min_pair_distance) - m_popts.ideal_pair_distance) << "\n";
  }

  const auto& best_pair =
      (best_aligns.size() == 1) ? best_aligns[0] : best_aligns[m_rr_idx++ % best_aligns.size()];

  if (brief_dbg) {
    std::cerr << "Saving to left:\n";
    dump_align(*best_pair.first, read_id, read_len);
    std::cerr << "Saving to right:\n";
    dump_align(*best_pair.second, rc_mate_id, mate_len);
  }
  save_align(*best_pair.first, read_id, read_len);
  save_align(*best_pair.second, rc_mate_id, mate_len);
}  // namespace variants

bool place_pair_cov::metric_better(const pair_align_metric& lhs,
                                   const pair_align_metric& rhs) const {
  int lhs_abs_dist = abs(lhs.dist_from_ideal);
  int rhs_abs_dist = abs(rhs.dist_from_ideal);

  if (lhs_abs_dist != rhs_abs_dist) {
    return lhs_abs_dist < rhs_abs_dist;
  }
  return false;
}

bool place_pair_cov::state_has_read(trace_state* st, uint32_t read_id, aoffset_t offset,
                                    int read_len) const {
  return st->a->read_coverage->get_read_ids_at(offset, read_len).contains(read_id);
}

void place_pair_cov::save_align(const align& aln, uint32_t read_id, int read_len) {
  for (const auto& part : aln.parts) {
    part.st->filtered_coverage.insert(part.offset, read_id, read_len);
  }
}

void place_pair_cov::dump_align(const align& aln, uint32_t read_id, int read_len) const {
  bool first = true;
  for (const auto& part : aln.parts) {
    const auto& a = *part.st->a;
    if (first) {
      first = false;
      std::cerr << " ref@" << a.left_offset << "\n";
    }
    std::cerr << "  AID " << a.assembly_id;
    if (a.matches_reference) {
      std::cerr << " REF len=" << a.right_offset - a.left_offset;
    } else {
      std::cerr << " reflen=" << a.right_offset - a.left_offset << " varlen=" << a.seq.size();
    }
    std::cerr << " rd@" << part.offset << "-" << part.offset + read_len << " ref@" << a.right_offset
              << "\n";
  }
}

bool place_pair_cov::propagate_read(
    uint32_t read_id, trace_state* st, aoffset_t offset, aoffset_t read_len,
    bool prop_left /* true if propagating left, false if propagating right */,
    std::vector<align>& align_out, align& so_far) {
  if (k_dbg) {
    std::cerr << "Attempting to propagate " << read_id << " to " << *st->a
              << " at offset=" << offset << " read_len= " << read_len
              << " dir=" << (prop_left ? "LEFT" : "RIGHT") << "\n";
  }
  aoffset_t seqlen = st->a->seq.size();
  aoffset_t end_offset = offset + read_len;

  CHECK_GT(end_offset, 0);
  CHECK_LT(offset, seqlen);
  DCHECK(state_has_read(st, read_id, offset, read_len));

  align_part new_part;
  new_part.st = st;
  new_part.offset = offset;
  so_far.parts.push_back(new_part);
  auto cleanup = [&]() {
    CHECK_EQ(so_far.parts.back(), new_part);
    so_far.parts.pop_back();
  };

  if (prop_left) {
    if (offset >= 0) {
      // Propagate complete
      align_out.push_back(so_far);
      cleanup();
      if (k_dbg) {
        std::cerr << "Finished propagating left\n";
      }
      return true;
    }
  } else {  // right
    if (end_offset <= seqlen) {
      // Propagate complete
      align_out.push_back(so_far);
      cleanup();
      if (k_dbg) {
        std::cerr << "Finished propagating right\n";
      }
      return true;
    }
  }

  const auto& edges = prop_left ? st->left_edges : st->right_edges;

  std::vector<trace_state*> next;

  for (trace_state* new_st : edges) {
    aoffset_t new_offset;
    if (prop_left) {
      aoffset_t new_seqlen = new_st->a->seq.size();
      new_offset = offset + new_seqlen;
    } else {
      new_offset = offset - seqlen;
    }
    if (k_dbg) {
      std::cerr << "Checking for read, new_offset= " << new_offset << ":" << *new_st->a << "\n";
    }
    if (state_has_read(new_st, read_id, new_offset, read_len)) {
      next.push_back(new_st);
    }
  }

  if (next.empty()) {
    cleanup();
    return false;
  }

  size_t rr_idx;
  if (prop_left) {
    rr_idx = st->left_edge_rr_idx++;
  } else {
    rr_idx = st->right_edge_rr_idx++;
  }

  size_t num_nexts = next.size();
  auto it = next.begin() + (rr_idx % num_nexts);

  for (size_t i = 0; i != num_nexts; ++i, ++it) {
    if (it == next.end()) {
      it = next.begin();
    }

    trace_state* new_st = *it;

    aoffset_t new_offset;
    if (prop_left) {
      aoffset_t new_seqlen = new_st->a->seq.size();
      new_offset = offset + new_seqlen;
    } else {
      new_offset = offset - seqlen;
    }

    if (propagate_read(read_id, new_st, new_offset, read_len, prop_left, align_out, so_far)) {
      cleanup();
      return true;
    }
  }
  cleanup();
  return false;
}

void place_pair_cov::save_filtered_reads() {
  for (auto& elem : m_asm_to_state) {
    auto* a = elem.first;
    auto& st = elem.second;

    a->pair_read_coverage.emplace(st->filtered_coverage.build_and_clear(a->seq.size()));
    st.reset();
  }
}

absl::btree_set<aoffset_t> place_pair_cov::testing_distances_between(aoffset_t left_offset,
                                                                     aoffset_t right_offset) const {
  auto right_it = m_dists.find(right_offset);
  if (right_it == m_dists.end()) {
    return {};
  }

  auto left_it = right_it->second.find(left_offset);
  if (left_it == right_it->second.end()) {
    return {};
  }

  const dist_set& dists = left_it->second;
  return absl::btree_set<aoffset_t>(dists.begin(), dists.end());
}

std::ostream& operator<<(std::ostream& os, const place_pair_cov::trace_state& st) {
  os << "Trace state for " << *st.a << " ";

  os << "left:";
  if (st.left_edges.empty()) {
    os << " (none)";
  } else {
    for (place_pair_cov::trace_state* other_st : st.left_edges) {
      if (other_st) {
        os << " " << other_st->a->assembly_id << "";
      } else {
        os << " (nullptr)";
      }
    }
  }
  os << "right:";
  if (st.right_edges.empty()) {
    os << " (none)";
  } else {
    for (place_pair_cov::trace_state* other_st : st.right_edges) {
      if (other_st) {
        os << " " << other_st->a->assembly_id << "";
      } else {
        os << " (nullptr)";
      }
    }
  }
  return os << "\n";
}

std::ostream& operator<<(std::ostream& os, const place_pair_cov::pair_align_metric& metric) {
  return os << "[dist=" << metric.dist_from_ideal << "]";
}

std::ostream& operator<<(std::ostream& os, const place_pair_cov::align_part& part) {
  return os << part.offset << "@" << *part.st->a;
}

std::ostream& operator<<(std::ostream& os, const place_pair_cov::read_info& ri) {
  if (ri.read_states.empty()) {
    os << "(no fwd states)";
  } else {
    os << "fwd states";
    for (const auto& rd : ri.read_states) {
      os << " " << rd.st->a->assembly_id << "@" << rd.st->a->right_offset << "+" << rd.offset;
    }
  }
  os << " ";

  if (ri.rc_mate_states.empty()) {
    os << "(no rc mate states)";
  } else {
    os << "rc mate states";
    for (const auto& rd : ri.rc_mate_states) {
      os << " " << rd.st->a->assembly_id << "@" << rd.st->a->left_offset << "+" << rd.offset;
    }
  }
  return os;
}

struct place_pair_cov::dist_set_formatter {
  void operator()(std::string* out, const dist_set_t& dists) const {
    out->append(StrJoin(dists, ","));
  }
  void operator()(std::string* out, const dist_set& dists) const {
    std::stringstream os;
    os << dists;
    out->append(os.str());
  }
};

std::string place_pair_cov::dump_dist_table(const dists_t& dt) const {
  std::stringstream os;
  os << absl::StrJoin(dt, " ",
                      absl::PairFormatter(absl::AlphaNumFormatter(), ":", dist_set_formatter()))
     << "\n";
  return os.str();
}

void place_pair_cov::dump_state(const std::string& where) {
  if (!k_dbg) {
    return;
  }
  std::cerr << "\nDumping placer state at '" << where << "', asms:\n";
  for (const auto& a : m_asm_to_state) {
    if (a.second) {
      std::cerr << *a.second << "\n";
    } else {
      std::cerr << *a.first << "\n";
    }
  }
  std::cerr << m_dists.size() << " dists:\n";
  for (const auto& dist : m_dists) {
    std::cerr << "right offset: " << dist.first << " dist table=" << dump_dist_table(dist.second)
              << "\n";
  }

  std::cerr << m_reads.size() << " reads:\n";
  for (const auto& rd : m_reads) {
    std::cerr << "read id=" << rd.first << ": " << rd.second << "\n";
  }
}

void place_pair_cov::on_assembly(assembly_ptr a) {
  auto new_state = make_unique<trace_state>(m_next_trace_id, a.get());
  ++m_next_trace_id;
  trace_state* st = new_state.get();

  bool did_insert = m_asm_to_state.insert(std::make_pair(a.get(), std::move(new_state))).second;
  CHECK(did_insert) << "Duplicate assembly " << *a << "?";

  save_reads(a.get(), st);

  m_block.emplace_back(std::move(a));
}

}  // namespace variants

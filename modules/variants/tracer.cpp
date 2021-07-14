#include "modules/variants/tracer.h"

#include <boost/range/adaptor/reversed.hpp>
#include <chrono>
#include <ctime>
#include <random>

#include "modules/bio_base/readmap.h"
#include "modules/io/parallel.h"

// Design considerations:
// Chromosome 1 is about 250 megabasepairs long.
// Regions that contain less than 10,000 N's in a row
// have a maximum length of 88 megabasepairs.

namespace variants {

const char tracer::k_tracer_name[] = "PUSH";

constexpr int k_asm_debug = 0;

static int tracer_debug_at(aoffset_t offset) {
  if (k_asm_debug) {
    return k_asm_debug;
  }
  if (offset_needs_trace(offset)) {
    return 3;
  }
  return 0;
}

#define NP_DEBUG(NP, EXPR)        \
  do {                            \
    if (m_options.debug_paths) {  \
      path* p = (NP).new_path;    \
      PATH_DEBUG(p, d.np = (NP)); \
      PATH_DEBUG(p, EXPR);        \
    }                             \
  } while (0)

#define PATH_DEBUG(P, EXPR)       \
  do {                            \
    if (m_options.debug_paths) {  \
      auto& d = m_path_debugs[P]; \
      (EXPR);                     \
    }                             \
  } while (0)

#define ADD_COST(NP, ACOST_TYPE, MULTIPLIER, SHORT_NAME)                                   \
  do {                                                                                     \
    /* Make sure max cost gets nowhere near maxint */                                      \
    CHECK_LT((NP).new_path->cost, (std::numeric_limits<acost_t>::max() / 2) / MULTIPLIER); \
    (NP).new_path->cost += m_options.ACOST_TYPE##_cost * MULTIPLIER;                       \
    ++m_stats.ACOST_TYPE##_cost;                                                           \
    NP_DEBUG(NP, d.costs[SHORT_NAME] += MULTIPLIER);                                       \
  } while (0)

class tracer::path_storage {
 public:
  path* alloc() {
    if (!m_storage || m_size == m_used) {
      CHECK_GT(m_size, 0);
      if (m_storage) {
        m_to_free.emplace_back(std::move(m_storage));
        m_size *= 2;
      }
      m_storage.reset(new path[m_size]);
      m_used = 0;
    }

    CHECK_LT(m_used, m_size);
    path* new_path(&m_storage[m_used]);
    ++m_used;

    *new_path = path();
    return new_path;
  }

  void reset() {
    m_used = 0;
    m_to_free.clear();
  }

  void migrate_to(path_storage& target) {
    for (auto& p : m_to_free) {
      target.m_to_free.emplace_back(std::move(p));
    }
    m_to_free.clear();
    if (m_storage) {
      target.m_to_free.emplace_back(std::move(m_storage));
    }
    m_storage.reset();
  }

 private:
  size_t m_size = 1 << 8;
  size_t m_used = 0;
  std::unique_ptr<path[]> m_storage;
  std::vector<std::unique_ptr<path[]>> m_to_free;
};

tracer::tracer(const assemble_options& options)
    : m_options(options),
      m_seqset(options.seqset),
      m_readmap(options.readmap),
      m_ref(options.ref),
      m_rmap(options.rmap) {
  m_trace_path_storage = make_unique<path_storage>();

  CHECK(m_seqset);
  CHECK(m_readmap);
  CHECK(m_ref);
  CHECK(m_rmap);
  CHECK(m_options.scaffold);

  m_ref_parts = &m_options.scaffold->extents();
  CHECK(!m_ref_parts->empty());
}

std::string tracer::path_to_string(const path* p) const {
  if (!p) {
    return "NOPATH";
  }

  std::stringstream result;
  result << "Entry=";
  if (p->range.begin() + 1 != p->range.end()) {
    result << "[" << p->range.begin() << "," << p->range.end() << ")";
  } else {
    result << p->range.begin();
  }
  result << ":" << int(p->range.size());
  result << " bases=" << p->seq.size() << ": " << p->seq;
  auto prev = p->prev;
  if (prev) {
    result << " prev: ";
    if (prev->range.begin() + 1 != prev->range.end()) {
      result << "[" << prev->range.begin() << "," << prev->range.end() << ")";
    } else {
      result << prev->range.begin();
    }
    result << ":" << int(prev->range.size());
  }
  return result.str();
}

std::string tracer::rejoin_to_string(const rejoin& r) const {
  std::stringstream result;
  result << "Rejoin " << r.p->cost << "+" << r.rejoin_cost;
  result << " at " << r.right_offset << "(" << r.right_anchor_len << "):";
  result << path_to_string(r.p);
  return result.str();
}

std::string tracer::next_path_to_string(const next_path& np) const {
  std::stringstream result;
  result << "[NP: min_overlap=" << int(np.new_path->min_overlap) << " a=" << np.ambiguous_bases
         << " br=" << np.branch_count_since_pair << " pairs: " << np.pairs_used
         << " path cost: " << np.new_path->cost << "]: path=" << path_to_string(np.new_path);
  return result.str();
}

struct tracer::next_path_comparer {
  // Compares next paths for the priority queue.  Returns true if rhs
  // is better than p.
  bool operator()(const next_path& p, const next_path& rhs) const {
    return p.new_path->cost > rhs.new_path->cost;
  }
};

// Compares rejoin candidates.  Returns true if rhs is better than p.
bool tracer::rejoin_comparer::operator()(const rejoin& p, const rejoin& rhs) const {
  return p.p->cost + p.rejoin_cost > rhs.p->cost + rhs.rejoin_cost;
}

struct tracer::path_debug_info {
  next_path np;
  boost::optional<rejoin> r;
  boost::optional<aoffset_t> ref_pos;
  seqset_range original_ref_range;
  std::vector<size_t> assembly_ids;
  std::vector<size_t> head_assembly_ids;
  std::map<std::string, acost_t> costs;
  bool explored = false;
  bool has_next = false;

  std::vector<std::string> tags;
  const path* prev_ref = nullptr;
};

assemble_stats tracer::assemble(assemble_pipeline_interface* output, progress_handler_t progress) {
  return assemble(0, std::numeric_limits<aoffset_t>::max(), output, progress);
}

assemble_stats tracer::assemble(aoffset_t start_offset, aoffset_t limit_offset,
                                assemble_pipeline_interface* output, progress_handler_t progress) {
  CHECK_LE(limit_offset, start_offset + aoffset_t(m_options.scaffold_split_size))
      << "start: " << start_offset << " limit: " << limit_offset;

  m_stats = assemble_stats();
  m_start_offset = start_offset;
  m_limit_offset = limit_offset;
  m_read_ahead_it = m_options.scaffold->begin();
  m_read_ahead_offset = m_read_ahead_it.offset();
  m_position_entry_index = 0;

  static std::mutex debug_mu;
  std::unique_lock<std::mutex> debug_l(debug_mu, std::defer_lock);
  if (k_asm_debug) {
    // Only single thread if we're outputting debug information.
    debug_l.lock();
    std::cout << "Starting with part at offset " << m_read_ahead_offset << "\n";
  }
  m_read_ahead_range = m_seqset->ctx_begin();

  skip_ahead_to(start_offset - m_options.max_pair_distance - m_options.seqset->max_read_len());

  size_t tot_bases = 0;
  for (const auto& part : *m_ref_parts) {
    if (part.offset >= limit_offset) {
      break;
    }
    if (part.offset + aoffset_t(part.sequence.size()) <= start_offset) {
      continue;
    }
    tot_bases += part.sequence.size();
  }
  size_t tot_bases_so_far = 0;
  for (const auto& part : *m_ref_parts) {
    if (part.offset < m_cur_offset) {
      throw(io_exception("Ref parts must be sorted"));
    }
    m_cur_offset = part.offset;
    if ((part.offset + aoffset_t(part.sequence.size())) <= start_offset) {
      continue;
    }

    if (part.offset > limit_offset) {
      break;
    }

    dna_slice to_process = part.sequence;
    if (part.offset + int(m_options.seqset->max_read_len()) < start_offset) {
      aoffset_t skip_len = start_offset - part.offset - int(m_options.seqset->max_read_len());
      CHECK_LT(skip_len, part.sequence.size()) << "Should have already skipped this one!";
      to_process = to_process.subseq(skip_len, aoffset_t(to_process.size()) - skip_len);
      m_cur_offset += skip_len;
      tot_bases_so_far += skip_len;
    }

    m_cur_left_anchor = dna_sequence();
    m_cur_ref_ambiguous_bases = 0;

    for (dna_base b : to_process) {
      advance_read_ahead_to(m_cur_offset + m_options.read_ahead_distance);
      advance_trail_behind_to(m_cur_offset - aoffset_t(m_options.max_pair_distance));

      m_cur_left_anchor.push_back(b);
      ++m_cur_ref_ambiguous_bases;
      if (m_cur_offset >= limit_offset) {
        return m_stats;
      }
      ++m_cur_offset;
      ++tot_bases_so_far;
      if ((tot_bases_so_far & 0xFF) == 0) {
        progress(tot_bases_so_far * 1.0 / tot_bases);
      }

      if (m_cur_offset < start_offset) {
        continue;
      }

      if (offset_needs_trace(m_cur_offset)) {
        std::cout.flush();
        std::cout << "Tracer got to offset " << m_cur_offset << "\n";
        std::cout.flush();
      }
      advance_position_entry_index();
      if (m_position_entry_index < m_position_entries.size() &&
          m_position_entries[m_position_entry_index].first == m_cur_offset) {
        m_cur_range = m_seqset->ctx_entry(m_position_entries[m_position_entry_index].second);

        if (is_read(m_cur_range)) {
          if (m_cur_left_anchor.size() < m_cur_range.size()) {
            continue;
          }
          m_cur_left_anchor = m_cur_left_anchor.subseq(
              m_cur_left_anchor.size() - m_cur_range.size(), m_cur_range.size());

          assemble_stats old_assemble_stats = m_stats;
          m_stats = assemble_stats();
          auto start_time = std::chrono::high_resolution_clock::now();
          trace();
          auto end_time = std::chrono::high_resolution_clock::now();
          double duration = std::chrono::duration<double>(end_time - start_time).count();
          static std::atomic<double> max_duration{0};
          double prev_max = max_duration.load();
          if (m_options.report_long_traces_func && duration > prev_max &&
              max_duration.compare_exchange_strong(prev_max, duration)) {
            m_options.report_long_traces_func(m_options.scaffold_name, duration, m_cur_offset,
                                              m_stats);
          }
          m_stats += old_assemble_stats;

          output_rejoins(output);
        }
      }

      if (offset_needs_trace(m_cur_offset)) {
        std::cout.flush();
        std::cout << "Tracer done tracing at offset " << m_cur_offset << "\n";
        std::cout.flush();
      }
    }
  }

  return m_stats;
}

void tracer::advance_position_entry_index() {
  while (m_position_entry_index != m_position_entries.size() &&
         m_position_entries[m_position_entry_index].first < m_cur_offset) {
    ++m_position_entry_index;
  }
}

void tracer::skip_ahead_to(aoffset_t offset) {
  if (offset <= m_read_ahead_offset) {
    return;
  }
  if (m_read_ahead_it == m_options.scaffold->end()) {
    return;
  }

  if (k_asm_debug) {
    std::cout << "Skipping from " << m_read_ahead_it.offset() << "(" << m_read_ahead_offset
              << ") to " << offset << "\n";
  }
  m_read_ahead_it.skip_to(offset, "tracer");
  m_read_ahead_offset = m_read_ahead_it.offset();
  m_read_ahead_range = m_seqset->ctx_begin().push_front_drop(m_read_ahead_it->complement());
  if (k_asm_debug) {
    std::cout << "Skipping ended up at " << m_read_ahead_it.offset() << "\n";
  }
}

void tracer::advance_read_ahead_to(aoffset_t offset) {
  while (m_read_ahead_offset < offset) {
    if (m_read_ahead_it == m_options.scaffold->end()) {
      m_read_ahead_offset = offset;
      return;
    }

    if (offset <= m_read_ahead_it.offset()) {
      return;
    }

    if (m_read_ahead_it.first_in_extent()) {
      m_read_ahead_range = m_seqset->ctx_begin();
      CHECK_GE(m_read_ahead_it.offset(), m_read_ahead_offset);
      m_read_ahead_offset = m_read_ahead_it.offset();
    } else {
      CHECK_EQ(m_read_ahead_offset, m_read_ahead_it.offset());
    }

    m_read_ahead_range = m_read_ahead_range.push_front_drop(m_read_ahead_it->complement());
    // We want the next offset after this base, before m_read_ahead_it
    // might possibly skip to the next extent.
    ++m_read_ahead_offset;
    ++m_read_ahead_it;

    if (is_read(m_read_ahead_range)) {
      uint64_t seqset_id = m_read_ahead_range.seqset_id();

      m_position_entries.push_back(std::make_pair(m_read_ahead_offset, seqset_id));

      bool skip_pairing = false;
      if (m_options.ignore_ambiguous_ref_pairs && m_rmap->get(seqset_id).match_count() > 1) {
        skip_pairing = true;
      }
      if (!skip_pairing) {
        auto read_ids = m_readmap->entry_to_index(seqset_id);
        for (auto read_id = read_ids.first; read_id != read_ids.second; ++read_id) {
          if (!m_readmap->has_mate(read_id)) {
            continue;
          }
          auto rc_mate_read_id = m_readmap->get_rev_comp(m_readmap->get_mate(read_id));
          m_rc_mate_read_positions.insert(std::make_pair(rc_mate_read_id, m_read_ahead_offset));
        }
      }
      m_entry_positions.insert(std::make_pair(seqset_id, m_read_ahead_offset));
    }

    if (m_read_ahead_it == m_options.scaffold->end()) {
      if (k_asm_debug) {
        std::cout << "No more parts\n";
      }
      return;
    }
  }
}

void tracer::advance_trail_behind_to(aoffset_t offset) {
  while (!m_position_entries.empty()) {
    aoffset_t old_offset = m_position_entries.front().first;
    uint64_t old_seqset_id = m_position_entries.front().second;

    if (old_offset >= offset) {
      return;
    }

    m_position_entries.pop_front();
    if (m_position_entry_index) {
      --m_position_entry_index;
    }

    bool skip_pairing = false;
    if (m_options.ignore_ambiguous_ref_pairs && m_rmap->get(old_seqset_id).match_count() > 1) {
      skip_pairing = true;
    }

    if (!skip_pairing) {
      auto read_ids = m_readmap->entry_to_index(old_seqset_id);
      for (auto read_id = read_ids.first; read_id != read_ids.second; ++read_id) {
        if (!m_readmap->has_mate(read_id)) {
          continue;
        }
        auto rc_mate_read_id = m_readmap->get_rev_comp(m_readmap->get_mate(read_id));
        auto eq_range = m_rc_mate_read_positions.equal_range(rc_mate_read_id);
        bool found = false;
        for (auto it = eq_range.first; it != eq_range.second; ++it) {
          if (it->second == old_offset) {
            m_rc_mate_read_positions.erase(it);
            found = true;
            break;
          }
        }
        CHECK(found);
      }
    }

    auto eq_range = m_entry_positions.equal_range(old_seqset_id);
    bool found = false;
    for (auto it = eq_range.first; it != eq_range.second; ++it) {
      if (it->second == old_offset) {
        m_entry_positions.erase(it);
        found = true;
        break;
      }
    }
    CHECK(found);
  }
}

bool tracer::has_seqset_id_in_range(uint64_t seqset_id, aoffset_t start_offset,
                                    aoffset_t limit_offset) const {
  auto eq_range = m_entry_positions.equal_range(seqset_id);
  for (auto it = eq_range.first; it != eq_range.second; ++it) {
    if (it->second >= start_offset && it->second < limit_offset) {
      return true;
    }
  }
  return false;
}

bool tracer::has_rc_mate_in_range(uint32_t read_id, aoffset_t start_offset,
                                  aoffset_t limit_offset) const {
  auto eq_range = m_rc_mate_read_positions.equal_range(read_id);
  for (auto it = eq_range.first; it != eq_range.second; ++it) {
    if (it->second >= start_offset && it->second < limit_offset) {
      return true;
    }
  }
  return false;
}

void tracer::trace() {
  if (!m_options.debug_paths) {
    m_trace_path_storage->reset();
    m_prev_path = nullptr;
  }
  m_next_paths.clear();
  m_rejoin_paths.clear();
  m_dead_end_rejoins.clear();
  m_trace_outputs = 0;

  if (m_cur_offset < m_start_offset || m_cur_offset >= m_limit_offset) {
    return;
  }

  bool ref_is_ambiguous = m_rmap->get(m_cur_range.seqset_id()).match_count() > 1;
  if (ref_is_ambiguous) {
    ++m_stats.ambiguous_ref_reads;

    if (!m_options.trace_ambiguous_ref) {
      return;
    }
  } else {
    ++m_stats.ref_reads;
    m_cur_ref_ambiguous_bases = 0;
  }

  next_path np;
  np.new_path = m_trace_path_storage->alloc();
  np.new_path->range = m_cur_range;
  np.new_path->seq = m_cur_left_anchor;
  np.ambiguous_bases = m_cur_ref_ambiguous_bases;
  np.new_path->min_overlap = m_cur_left_anchor.size();

  NP_DEBUG(np, d.prev_ref = m_prev_path);
  NP_DEBUG(np, d.original_ref_range = m_cur_range);
  NP_DEBUG(np, d.ref_pos.emplace(m_cur_offset));
  m_prev_path = np.new_path;

  push_next_path(std::move(np));

  trace_next_paths();
}

bool tracer::is_read(const seqset_range& r) const {
  if (!r.is_seqset_entry()) {
    return false;
  }

  return m_readmap->get_bit(r.seqset_id());
}

void tracer::expand_next_path(next_path np) {
  path* new_path = np.new_path;

  unsigned valid_count = 0;
  dna_base_array<seqset_range> pushed;
  bool pushed_any = false;

  NP_DEBUG(np, d.explored = true);

  if (tracer_debug_at(m_cur_offset) > 3) {
    std::cout << "Expanding path: " << next_path_to_string(np) << "\n";
  }

  // Add the next base along this path
  std::multiset<unsigned> push_lengths;
  for (;;) {
    push_lengths.clear();
    valid_count = 0;
    dna_base next_base{};
    if (tracer_debug_at(m_cur_offset) > 3) {
      std::cout << "Requiring " << m_options.min_overlap + np.pushed_since_read
                << " bases because of " << np.pushed_since_read
                << " pushed so far, pushing onto: " << np.new_path->range.sequence() << "\n";
    }
    for (dna_base b : dna_bases()) {
      pushed[b] = np.new_path->range.push_front_drop(b.complement(),
                                                     m_options.min_overlap + np.pushed_since_read);
      if (!pushed[b].valid()) {
        continue;
      }
      if (tracer_debug_at(m_cur_offset) > 3) {
        std::cout << "Pushing " << b << " results in: " << pushed[b].sequence() << "\n";
      }
      if (tracer_debug_at(m_cur_offset) > 3) {
        if (is_read(pushed[b])) {
          std::cout << "Found entry at base " << b << ": " << pushed[b].seqset_id() << ": "
                    << pushed[b].sequence() << "\n";
        }
      }
      if (pushed[b].is_seqset_entry()) {
        if (np.loop_check_seqset_id == pushed[b].seqset_id()) {
          ++m_stats.loops;
          pushed[b] = seqset_range();
          if (tracer_debug_at(m_cur_offset) > 3) {
            std::cout << "Skipping loop: " << pushed[b].seqset_id() << "\n";
          }
          continue;
        }
      }
      valid_count++;
      next_base = b;
      push_lengths.insert(pushed[b].size());
    }
    if (valid_count != 1) {
      ++m_ambiguous_search_step_count;
      break;
    }
    if (pushed_any && is_read(new_path->range)) {
      if (tracer_debug_at(m_cur_offset) > 3) {
        std::cout << "Found read; done expanding this trace\n";
      }
      break;
    }

    const auto& next_pushed = pushed[next_base];
    CHECK(next_pushed.valid());
    new_path->range = next_pushed;
    if (tracer_debug_at(m_cur_offset) > 3) {
      std::cout << next_path_to_string(np) << " has single path, range: " << next_pushed.sequence()
                << " next_base: " << next_base << "\n";
    }
    add_base_to_next_path(np, next_base);
    pushed_any = true;

    if (np.pushed_since_pair > m_options.max_bases_between_pairs) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "DISCARD PATH: too far without pair; " << np.pushed_since_pair << " > "
                  << m_options.max_bases_between_pairs << ": " << next_path_to_string(np) << "\n";
      }
      NP_DEBUG(np, d.tags.push_back("too-far-without-pair"));
      ++m_stats.too_far_without_pair;
      return;
    }

    if (add_rejoins(np)) {
      if (tracer_debug_at(m_cur_offset) > 3) {
        std::cout << "Done extending; Rejoins added.\n";
      }
      NP_DEBUG(np, d.tags.push_back("added-rejoins-1"));
      return;
    }
  }

  if (tracer_debug_at(m_cur_offset) > 3) {
    std::cout << push_lengths.size() << " options from here\n";
  }

  if (new_path->prev) {
    add_dead_end_rejoin(np);
  }

  if (push_lengths.empty()) {
    if (tracer_debug_at(m_cur_offset)) {
      std::cout << "DISCARD PATH: no paths forward: " << next_path_to_string(np) << "\n";
    }

    NP_DEBUG(np, d.tags.push_back("no-paths-forward"));
    return;
  }

  auto it = push_lengths.rbegin();
  unsigned unambiguous_len = *it;
  ++it;
  if (it != push_lengths.rend() && *it == unambiguous_len) {
    // Two with the same length; all paths are ambiguous.
    ++unambiguous_len;
  }

  for (dna_base b : dna_bases()) {
    if (pushed[b].valid()) {
      next_path new_np = np;
      CHECK_EQ(new_path, new_np.new_path);
      new_np.new_path = m_trace_path_storage->alloc();
      new_np.new_path->prev = new_path;
      new_np.new_path->min_overlap = new_path->min_overlap;
      new_np.new_path->cost = new_path->cost;
      new_np.new_path->range = pushed[b];
      add_base_to_next_path(new_np, b);
      if (pushed[b].size() < unambiguous_len) {
        ADD_COST(new_np, ambiguous_branch, 1, "?br");
        ++new_np.branch_count_since_pair;
        if (new_np.branch_count_since_pair > m_options.max_branches_between_pairs) {
          ++m_stats.exceeded_branch_limit;
          NP_DEBUG(new_np, d.tags.push_back("max-branch-count-since-pair"));
          continue;
        }
      }
      NP_DEBUG(new_np, d.explored = false);
      if (add_rejoins(new_np)) {
        NP_DEBUG(new_np, d.tags.push_back("added-rejoins-2"));
      } else {
        push_next_path(std::move(new_np));
      }
    }
  }
}

bool tracer::add_rejoins(const next_path& np) {
  if (!is_read(np.new_path->range)) {
    return false;
  }
  uint64_t seqset_id = np.new_path->range.seqset_id();

  auto eq_range = m_entry_positions.equal_range(seqset_id);
  aoffset_t best_rejoin_distance = std::numeric_limits<aoffset_t>::max();
  aoffset_t best_rejoin_offset = -1;
  for (auto it = eq_range.first; it != eq_range.second; ++it) {
    if (it->second <= m_cur_offset) {
      continue;
    }
    aoffset_t ideal_pos = m_cur_offset + np.path_bases;
    aoffset_t candidate_pos = it->second;
    aoffset_t distance = llabs(ideal_pos - candidate_pos);
    if (distance < best_rejoin_distance) {
      best_rejoin_distance = distance;
      best_rejoin_offset = candidate_pos;
    }
  }

  if (best_rejoin_distance != std::numeric_limits<aoffset_t>::max()) {
    if (np.ambiguous_bases && best_rejoin_distance > aoffset_t(np.ambiguous_bases)) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "Extending rejoin; " << best_rejoin_distance << " > " << np.ambiguous_bases
                  << ": " << next_path_to_string(np) << " offset: " << best_rejoin_offset
                  << " cur: " << m_cur_offset << " ideal: " << (m_cur_offset + np.path_bases)
                  << " diff: " << aoffset_t(m_cur_offset + np.path_bases) - best_rejoin_offset
                  << "\n";
      }
      NP_DEBUG(np, d.tags.push_back("ext-ambig-rj"));
      ++m_stats.extend_ambiguous_rejoin;
      return false;
    }

    rejoin new_rejoin;
    new_rejoin.p = np.new_path;
    new_rejoin.rejoin_cost = m_options.rejoin_local_cost;
    new_rejoin.right_offset = best_rejoin_offset;
    new_rejoin.right_anchor_len = np.new_path->range.size();
    new_rejoin.rejoin_cost += best_rejoin_distance * m_options.size_change_cost;
    new_rejoin.rejoin_cost -= m_options.traverse_ref_cost;
    push_rejoin(std::move(new_rejoin));

    return true;
  }

  return false;
}

void tracer::add_dead_end_rejoin(const next_path& np) {
  if (np.num_reads <= 1) {
    return;
  }
  if (!m_options.trace_dead_ends) {
    return;
  }
  m_dead_end_rejoins.insert(np.new_path);
  if (np.new_path->prev && np.new_path->prev->min_overlap == np.new_path->min_overlap) {
    m_dead_end_rejoins.erase(np.new_path->prev);
  }
}

void tracer::add_pairs_to_next_path(next_path& np) {
  if (!is_read(np.new_path->range)) {
    return;
  }

  uint64_t seqset_id = np.new_path->range.seqset_id();

  ++np.num_reads;

  if (np.reads_until_loop_check) {
    --np.reads_until_loop_check;
  } else {
    np.loop_check_seqset_id = seqset_id;
    np.reads_until_loop_check = np.num_reads / 2;
  }

  auto reads = m_readmap->entry_to_index(seqset_id);
  if (reads.first != reads.second) {
    np.pushed_since_read = 0;

    auto ref_count = m_rmap->get(seqset_id).match_count();
    if (ref_count > 1) {
      NP_DEBUG(np, d.tags.push_back("ambig-ref"));
      np.ambiguous_bases = std::max<unsigned>(np.ambiguous_bases, np.new_path->range.size());
      if (!m_options.trace_ambiguous_ref) {
        np.new_path->cost += m_options.max_cost;
        m_stats.prune_ambiguous_ref++;
      }
    } else if (ref_count == 1) {
      ADD_COST(np, traverse_ref, 1, "xactref");
    }
  }

  bool matched_pair = false;
  unsigned added_pairs = 0;
  unsigned added_reads = 0;
  for (auto read_id = reads.first; read_id != reads.second; ++read_id) {
    if (!m_readmap->has_mate(read_id)) {
      continue;
    }

    if (added_reads < m_options.max_pairs_per_read) {
      np.new_path->seen_read_ids.push_back(read_id);
      ++added_reads;
    }

    if (m_readmap->get_is_forward(read_id) == m_options.forward_pairs_face_inward) {
      if (matched_pair) {
        // Already matched a pair; no need to do it again
        continue;
      }

      // Pair to the left in reference:
      aoffset_t ambiguous_in_ref = 0;
      if (np.ambiguous_bases > np.path_bases) {
        ambiguous_in_ref = np.ambiguous_bases - np.path_bases;
      }
      if (has_rc_mate_in_range(read_id, m_cur_offset + np.path_bases - m_options.max_pair_distance,
                               m_cur_offset - std::max<aoffset_t>(0, ambiguous_in_ref))) {
        matched_pair = true;
        continue;
      }

      if (path_has_read_in_range(np.new_path, read_id, np.ambiguous_bases,
                                 m_options.max_pair_distance)) {
        matched_pair = true;
        continue;
      }
    } else {
      if (added_pairs < m_options.max_pairs_per_read) {
        path::seen_pair s;
        s.read_id = m_readmap->get_rev_comp(m_readmap->get_mate(read_id));
        s.offset = np.new_path->seq.size();
        np.new_path->seen_pairs.push_back(s);
        ++added_pairs;
        ++m_stats.found_pairs;
        if (added_pairs >= m_options.max_pairs_per_read) {
          ++m_stats.too_many_pairs;
        }
      }

      if (matched_pair) {
        continue;
      }

      // Pair to the right in reference:
      if (has_rc_mate_in_range(read_id, m_cur_offset,
                               m_cur_offset + np.path_bases + m_options.max_pair_distance)) {
        matched_pair = true;
        continue;
      }
    }
  }
  if (matched_pair) {
    np.pushed_since_pair = 0;
    np.branch_count_since_pair = 0;
    np.ambiguous_bases = 0;
    ++np.pairs_used;
    ADD_COST(np, pairs_used, 1, "pu");
    ++m_stats.matched_pairs;
  }
}

void tracer::add_base_to_next_path(next_path& np, dna_base b) {
  np.new_path->seq.push_back(b);
  ++np.pushed_since_read;
  ++np.path_bases;
  if (np.ambiguous_bases) {
    ++np.ambiguous_bases;
  }
  ADD_COST(np, base, 1, "b");
  CHECK_GT(np.new_path->range.size(), np.pushed_since_read) << next_path_to_string(np);
  CHECK_LT(np.pushed_since_read, np.new_path->range.size()) << next_path_to_string(np);
  unsigned overlap = np.new_path->range.size() - np.pushed_since_read;
  if (overlap < np.new_path->min_overlap) {
    CHECK_LT(np.new_path->min_overlap, std::numeric_limits<unsigned>::max());
    auto decrease_amt = np.new_path->min_overlap - overlap;
    ADD_COST(np, decrease_overlap, decrease_amt, "ol");
    np.new_path->min_overlap = overlap;
    // Give a free pair distance increase with overlap decrese; we
    // expect every other pair to give a confirmation, so multiply by
    // 2:
    np.max_between_pairs += decrease_amt * 2;
  }

  ++np.pushed_since_pair;
  if (np.pushed_since_pair > np.max_between_pairs) {
    ++np.max_between_pairs;
    ADD_COST(np, increase_max_between_pair, 1, "mbp");
    CHECK_EQ(np.max_between_pairs, np.pushed_since_pair)
        << "We should only increase maximum pair distance one base at a time";
  }

  add_pairs_to_next_path(np);
}

void tracer::push_next_path(next_path np) {
  if (np.new_path->cost > m_options.max_cost) {
    if (tracer_debug_at(m_cur_offset) > 1) {
      std::cout << "Max cost exceeded pushing next path " << next_path_to_string(np) << "\n";
    }
    NP_DEBUG(np, d.tags.push_back("exceeded-max-cost"));
    ++m_stats.max_branch_cost;
    return;
  }
  if (!m_rejoin_paths.empty() && m_rejoin_paths.size() >= m_options.max_rejoins) {
    auto last_rejoin = m_rejoin_paths.begin();
    unsigned worst_rejoin_cost = last_rejoin->p->cost + last_rejoin->rejoin_cost;
    if (np.new_path->cost > worst_rejoin_cost) {
      NP_DEBUG(np, d.tags.push_back("subopt-prune-2"));
      ++m_stats.suboptimal_path_prune;
      return;
    }
  }
  m_next_paths.emplace_back(std::move(np));
  std::push_heap(m_next_paths.begin(), m_next_paths.end(), next_path_comparer());
  if (m_next_paths.size() > m_options.max_next_paths) {
    NP_DEBUG(m_next_paths.back(), d.tags.push_back("queue-too-big"));
    m_next_paths.pop_back();
    ++m_stats.next_paths_too_big;
  }
}

tracer::next_path tracer::pop_next_path() {
  CHECK(!m_next_paths.empty());
  std::pop_heap(m_next_paths.begin(), m_next_paths.end(), next_path_comparer());
  next_path out = std::move(m_next_paths.back());
  m_next_paths.pop_back();
  return out;
}

void tracer::push_rejoin(rejoin r) {
  PATH_DEBUG(r.p, d.r.emplace(r));
  if (r.p->cost + r.rejoin_cost > m_options.max_cost) {
    if (tracer_debug_at(m_cur_offset) > 1) {
      std::cout << "Max cost exceeded pushing rejoin " << rejoin_to_string(r) << "\n";
    }
    ++m_stats.max_branch_cost;
    return;
  }
  if (tracer_debug_at(m_cur_offset) > 1) {
    std::cout << "Pushing rejoin: " << rejoin_to_string(r) << "\n";
  }
  CHECK(r.p);
  CHECK(is_read(r.p->range));
  m_rejoin_paths.emplace(std::move(r));
  if (m_rejoin_paths.size() > m_options.max_rejoins) {
    m_rejoin_paths.erase(m_rejoin_paths.begin());
  }
}

void tracer::trace_next_paths() {
  m_ambiguous_search_step_count = 0;
  m_search_step_count = 0;
  m_trace_outputs = 0;

  if (tracer_debug_at(m_cur_offset)) {
    std::cout << "Tracing " << m_next_paths.size() << " next paths starting at offset "
              << m_cur_offset << ", read ahead offset = " << m_read_ahead_offset
              << ", ambiguous bases = " << m_cur_ref_ambiguous_bases << "\n";
  }

  while (!m_next_paths.empty()) {
    next_path np = pop_next_path();

    // Filtering based on cost is free.
    if (!m_rejoin_paths.empty() && m_rejoin_paths.size() >= m_options.max_rejoins &&
        np.new_path->cost >
            (m_rejoin_paths.begin()->rejoin_cost + m_rejoin_paths.begin()->p->cost)) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "DISCARD PATH: Suboptimal path prune: cost " << np.new_path->cost
                  << " is worse than " << m_rejoin_paths.begin()->rejoin_cost << " + "
                  << m_rejoin_paths.begin()->p->cost << ": " << next_path_to_string(np) << "\n";
      }
      NP_DEBUG(np, d.tags.push_back("subo-drop"));
      ++m_stats.suboptimal_path_prune;
      continue;
    }

    if (m_search_step_count >
        (m_options.max_search_steps_per_read * np.num_reads + m_options.initial_search_steps)) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "DISCARD PATH: Ran out of search steps per read: " << next_path_to_string(np)
                  << "\n";
      }
      NP_DEBUG(np, d.tags.push_back("too-slow"));
      ++m_stats.search_not_fast_enough;
      continue;
    }

    if (np.ambiguous_bases > m_options.max_ambiguous_bases) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "DISCARD PATH: Too many ambiguous bases: " << next_path_to_string(np) << "\n";
      }
      NP_DEBUG(np, d.tags.push_back("too-much-ambig"));
      ++m_stats.too_many_ambiguous_bases;
      continue;
    }
    ++m_stats.step_count;
    ++m_search_step_count;

    if (m_search_step_count > m_options.max_search_steps) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "DISCARD PATH: Ran out of search steps: " << next_path_to_string(np) << "\n";
      }
      NP_DEBUG(np, d.tags.push_back("too-many-steps"));
      ++m_stats.too_many_steps;
      return;
    }
    if (m_ambiguous_search_step_count > m_options.max_ambiguous_search_steps) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "DISCARD PATH: Ran out of ambiguous search steps: " << next_path_to_string(np)
                  << "\n";
      }
      NP_DEBUG(np, d.tags.push_back("too-many-ambig-steps"));
      ++m_stats.too_many_ambiguous_steps;
      return;
    }

    if (tracer_debug_at(m_cur_offset)) {
      if (tracer_debug_at(m_cur_offset) > 1 || !(m_search_step_count % 1000)) {
        std::cout << "Step " << m_search_step_count << " (" << m_next_paths.size()
                  << " left): " << next_path_to_string(np) << "\n";
      }
    }

    expand_next_path(std::move(np));
  }
  m_stats.max_ambiguous_step_count.add(m_ambiguous_search_step_count);
}

void tracer::output_rejoins(assemble_pipeline_interface* output) {
  unsigned output_count = 0;

  if (tracer_debug_at(m_cur_offset) && !m_rejoin_paths.empty()) {
    std::cout << "Outputting " << m_rejoin_paths.size() << "/" << m_options.max_rejoins
              << " rejoins\n";
  }

  while (output_count < m_options.max_rejoins && !m_rejoin_paths.empty()) {
    auto last_it = m_rejoin_paths.end();
    --last_it;
    rejoin r = *last_it;
    m_rejoin_paths.erase(last_it);

    if (tracer_debug_at(m_cur_offset)) {
      std::cout << "Rejoin candidate: " << rejoin_to_string(r) << "\n";
    }
    if (r.p->part_of_assembly) {
      if (tracer_debug_at(m_cur_offset)) {
        std::cout << "Skipping because already output\n";
      }
      continue;
    }

    // Exhausted other options; output this dead end we found.
    output_assembly(r, output);
    ++output_count;
  }

  if (tracer_debug_at(m_cur_offset) && !m_rejoin_paths.empty()) {
    std::cout << "DISCARD " << m_rejoin_paths.size() << " rejoins due to max_rejoins\n";
    for (const auto& r : m_rejoin_paths) {
      std::cout << "DISCARD: " << rejoin_to_string(r) << "\n";
    }
  }
  if (output_count < m_options.max_rejoins) {
    output_dead_ends(m_options.max_rejoins - output_count, output);
  }
}

struct tracer::dead_end_sorter {
  bool operator()(const path* a, const path* b) const {
    // Must be a total ordering so that we don't have non-determinism.
    if (a->min_overlap != b->min_overlap) {
      return a->min_overlap > b->min_overlap;
    }
    if (a->cost != b->cost) {
      return a->cost < b->cost;
    }
    if (a->range != b->range) {
      if (a->range.size() != b->range.size()) {
        return a->range.size() > b->range.size();
      }
      return a->range.begin() < b->range.begin();
    }
    if (a->part_of_assembly != b->part_of_assembly) {
      return b->part_of_assembly;
    }
    if (a->seen_pairs != b->seen_pairs) {
      return a->seen_pairs < b->seen_pairs;
    }
    if (a->seen_read_ids != b->seen_read_ids) {
      return a->seen_read_ids < b->seen_read_ids;
    }
    if (a->seq != b->seq) {
      return a->seq < b->seq;
    }
    // Everything else is the same; compare the previous path.
    if (a->prev != b->prev) {
      if (bool(a->prev) != bool(b->prev)) {
        return bool(a->prev) < bool(b->prev);
      }

      return (*this)(a->prev, b->prev);
    }
    return false;
  }
};

void tracer::output_dead_ends(unsigned max_to_output, assemble_pipeline_interface* output) {
  unsigned output_count = 0;

  std::vector<const path*> dead_ends;
  dead_ends.insert(dead_ends.end(), m_dead_end_rejoins.begin(), m_dead_end_rejoins.end());
  m_dead_end_rejoins.clear();

  std::sort(dead_ends.begin(), dead_ends.end(), dead_end_sorter());

  if (tracer_debug_at(m_cur_offset)) {
    std::cout << "Considering " << dead_ends.size() << " dead end rejoins\n";
  }

  for (const path* p : dead_ends) {
    if (output_count >= max_to_output) {
      return;
    }

    if (p->part_of_assembly) {
      continue;
    }

    rejoin r;
    r.p = p;
    r.rejoin_cost = m_options.dead_end_cost + p->cost;
    r.anchor_drop = true;

    if (tracer_debug_at(m_cur_offset)) {
      std::cout << "Outputting dead end rejoin, cost " << p->cost << ": " << path_to_string(p)
                << " from " << m_cur_offset << "\n";
    }

    output_assembly(r, output);
    ++output_count;
  }
  if (tracer_debug_at(m_cur_offset)) {
    std::cout << "Done considering dead end rejoins\n";
  }
}

bool tracer::path_has_read_in_range(const path* start_path, uint32_t read_id, aoffset_t start,
                                    aoffset_t limit) const {
  const path* cur = start_path;
  CHECK(cur);

  aoffset_t cur_distance = 0;
  const path* cur_path = start_path;
  while (cur_path) {
    aoffset_t last_offset = start_path->seq.size();
    for (auto it = start_path->seen_pairs.rbegin(); it != start_path->seen_pairs.rend(); ++it) {
      path::seen_pair s = *it;
      CHECK_LE(s.offset, last_offset);
      cur_distance += (last_offset - s.offset);
      if (cur_distance >= limit) {
        return false;
      }
      if (cur_distance < start) {
        continue;
      }
      if (s.read_id == read_id && cur_distance > 0) {
        return true;
      }
    }
    cur_distance += last_offset;
    if (cur_distance >= limit) {
      return false;
    }

    cur_path = cur_path->prev;
  }

  return false;
}

void tracer::output_assembly(const rejoin& r, assemble_pipeline_interface* output) const {
  assembly_ptr out = make_unique<assembly>();
  out->tags.insert(k_tracer_name);
  out->assembly_id = allocate_assembly_id();
  PATH_DEBUG(r.p, d.head_assembly_ids.push_back(out->assembly_id));
  if (tracer_debug_at(m_cur_offset)) {
    std::cout << "Constructing assembly from " << m_cur_offset << "(-" << m_cur_left_anchor.size()
              << ")"
              << " to " << r.right_offset << "(" << r.right_anchor_len << "):";
  }
  const path* cur = r.p;
  CHECK(cur);

  dna_sequence seq;
  while (cur) {
    if (tracer_debug_at(m_cur_offset)) {
      std::cout << " cur->seq";
    }
    PATH_DEBUG(cur, d.assembly_ids.push_back(out->assembly_id));
    seq += cur->seq.rev_comp();
    cur->part_of_assembly = true;
    for (uint32_t read_id : cur->seen_read_ids) {
      out->rc_read_ids.insert(read_id);
    }
    // Don't let other assemblies reuse these for pair matching.
    cur->seen_read_ids.clear();
    cur = cur->prev;
  }

  if (tracer_debug_at(m_cur_offset)) {
    std::cout << "\n";
  }

  out->left_offset = m_cur_offset - aoffset_t(m_cur_left_anchor.size());
  out->left_anchor_len = m_cur_left_anchor.size();
  if (r.anchor_drop) {
    out->right_offset = m_cur_offset + seq.size() * m_options.anchor_drop_size_multiplier;
    out->right_anchor_len = 0;
    out->score += m_options.anchor_drop_score;
  } else {
    out->right_offset = r.right_offset;
    out->right_anchor_len = r.right_anchor_len;
    CHECK_LE(out->right_offset, m_options.scaffold->end_pos());
  }
  out->min_overlap = r.p->min_overlap;
  out->trace_steps = m_search_step_count;
  out->left_anchor_ambiguous_bases = m_cur_ref_ambiguous_bases;
  out->seq = seq.rev_comp();
  if (m_options.calculate_coverage) {
    out->coverage = m_readmap->approx_coverage(out->seq);
  }
  out->matches_reference =
      (out->left_anchor_len + out->right_anchor_len) >= (out->right_offset - out->left_offset) &&
      (out->right_offset - out->left_offset) == aoffset_t(out->seq.size());
  if (out->matches_reference) {
    out->left_anchor_len = out->right_anchor_len = 0;
  }
  CHECK_LE(out->left_offset, m_options.scaffold->end_pos());

  aoffset_t max_tot_anchor =
      std::min<aoffset_t>(out->seq.size(), out->right_offset - out->left_offset);
  if ((out->right_anchor_len + out->left_anchor_len) > max_tot_anchor) {
    // TODO(nils): Research this case and make sure it's well tested everywhere
    CHECK_GT(out->right_anchor_len, 0);
    out->right_anchor_len = std::min<aoffset_t>(out->right_anchor_len, max_tot_anchor / 2);
    out->left_anchor_len = std::min<aoffset_t>(out->left_anchor_len, max_tot_anchor / 2);
  }

  m_stats.max_assembly_len.add(out->seq.size());

  if (r.anchor_drop && m_options.report_anchor_drop_func) {
    m_options.report_anchor_drop_func(*out, false /* anchored on left */);
  }

  out->tags.insert("PUSH");

  if (assembly_needs_trace(*out)) {
    std::cout << "OUT: tracer " << this << " produced " << (void*)out.get() << ": " << *out << "\n";
  }
  output->add(std::move(out));
  ++m_stats.output_count;
}

tracer::~tracer() {
  if (m_options.debug_paths) {
    std::stringstream dot_contents;
    output_path_debug_dot(dot_contents);
    m_options.debug_paths(dot_contents.str());
  }
}

void tracer::output_path_debug_dot(std::ostream& os) const {
  std::vector<const path*> fills;
  for (const auto& pd : m_path_debugs) {
    const auto& p = pd.first;
    fills.push_back(p);
  }
  std::map<aoffset_t, std::string> ref_nodes;

  while (!fills.empty()) {
    auto p = fills.back();
    fills.pop_back();

    auto& d = m_path_debugs[p];
    if (d.ref_pos) {
      ref_nodes[*d.ref_pos] = printstring("P%p", p);
    }

    if (!p->prev) {
      continue;
    }
    if (!m_path_debugs.count(p->prev)) {
      fills.push_back(p->prev);
    }

    m_path_debugs[p->prev].has_next = true;
  }

  os << "digraph G {\n"
     << "  mode=\"hier\";"
     << "  ranksep=.2;\n"
     << "  node [shape=record, width=.1, height=.1];\n";

  for (const auto& pd : m_path_debugs) {
    const auto& p = pd.first;
    const auto& d = pd.second;

    std::string p_id = printstring("P%p", p);

    const path_debug_info* prev_d = nullptr;
    std::string prev_id;
    const path* prev = nullptr;
    if (p->prev) {
      prev = p->prev;
      prev_id = printstring("P%p", p->prev);
      prev_d = &m_path_debugs[p->prev];
    } else if (d.prev_ref) {
      prev = d.prev_ref;
      prev_id = printstring("P%p", d.prev_ref);
      prev_d = &m_path_debugs[d.prev_ref];
    }

    if (prev_d) {
      //      if (/* !d.explored &&  */ !d.ref_pos && prev_d->ref_pos &&
      //      !d.has_next) {
      // auto orig_ref_pos = *prev_d->ref_pos;
      // orig_ref_pos += p->prev->seq.size();
      // dna_slice ref_slice = m_rmap->get_ref_slice(orig_ref_pos);
      // if (ref_slice.size() > 100) {
      //   ref_slice = ref_slice.subseq(0, 100);
      // }
      // std::cout << "Comparing " << ref_slice << " to " << p->seq << "\n";
      // if (p->seq.compare_to(ref_slice) ==
      // dna_compare_result::FIRST_IS_PREFIX) {
      // Skip this one; it's probably just a node  to trace to the next ref
      // entry.
      //        continue;
      //        }
      //      }*/

      if (d.ref_pos) {
        os << "  " << prev_id << ":s -> " << p_id << ":n [weight=50, color=\"green\"]";
      } else {
        if (p->part_of_assembly) {
          if (prev_d->ref_pos) {
            os << "  " << prev_id << ":e -> " << p_id << ":w [weight=2, color=\"red\"]";
          } else {
            os << "  " << prev_id << ":s -> " << p_id << ":n [weight=4, color=\"red\"]";
          }
        } else {
          if (prev_d->ref_pos) {
            os << "  " << prev_id << ":w -> " << p_id;
          } else {
            os << "  " << prev_id << ":s -> " << p_id << ":n";
          }
        }
      }
      os << ";\n";
    }
    if (d.r && (!prev_d || !d.ref_pos)) {
      if (ref_nodes.count(d.r->right_offset)) {
        auto ref_id = ref_nodes[d.r->right_offset];
        os << "  " << p_id << ":rejoin_local:s -> " << ref_id << ":n [weight=1];\n";
      }
    }

    os << "  " << p_id << " [";

    if (d.ref_pos) {
      if (d.explored) {
        os << "color=\"green\", ";
      } else {
        os << "color=\"blue\", ";
      }
    } else if (p->part_of_assembly) {
      os << "color=\"red\", ";
    } else if (!d.explored) {
      os << "color=\"yellow\", ";
    }
    os << "label = \"{<seq> ";
    dna_slice show_seq = p->seq;
    if (d.ref_pos && d.prev_ref && d.original_ref_range.valid()) {
      aoffset_t seq_tail = aoffset_t(show_seq.size()) - aoffset_t(d.original_ref_range.size());
      CHECK_GE(seq_tail, 0);
      if (seq_tail && seq_tail < aoffset_t(show_seq.size())) {
        show_seq = show_seq.subseq(show_seq.size() - seq_tail, seq_tail);
      }
    }

    std::string show_seq_str;
    if (prev) {
      seqset_range r = d.original_ref_range.valid() ? d.original_ref_range : prev->range;
      if (is_read(r)) {
        show_seq_str += "!";
      }
      for (dna_base b : show_seq) {
        auto new_r = r.push_front_drop(b.complement());
        int dropped = int(r.size() + 1) - int(new_r.size());
        r = new_r;
        if (dropped) {
          show_seq_str += "(-" + std::to_string(dropped) + ")";
        }
        show_seq_str += char(b);
        if (is_read(new_r)) {
          show_seq_str += "!";
        }
      }
    } else {
      show_seq_str = "?" + show_seq.as_string();
    }

    if (show_seq_str.size() > 250) {
      os << show_seq_str.substr(0, 100) << "..." << show_seq_str.substr(show_seq_str.size() - 100);
    } else {
      os << show_seq_str;
    }
    os << " |{";

    std::string label_flags = "";
    aoffset_t seq_tail = 0;

    if (d.ref_pos) {
      label_flags += "<ref> Ref";
      if (seq_tail) {
        label_flags += printstring("%+d", seq_tail);
      }
      label_flags += ":" + std::to_string(*d.ref_pos);

      label_flags += " |";
    } else {
      aoffset_t path_bases = 0;
      const path* t = p;
      while (t && t->prev) {
        path_bases += t->seq.size();
        t = t->prev;
      }
      CHECK(m_path_debugs[t].ref_pos);
      aoffset_t rp = *m_path_debugs[t].ref_pos;
      label_flags += printstring("Ref+%d~%d |", path_bases, rp + path_bases);
    }

    if (p->part_of_assembly) {
      label_flags += " Asm |";
    }

    if (d.r) {
      if (d.r->right_anchor_len) {
        label_flags += " <rejoin_local> Rejoin at " + std::to_string(d.r->right_offset) + "(" +
                       std::to_string(d.r->right_anchor_len) +
                       ") cost=" + std::to_string(d.r->p->cost) +
                       " rjcost=" + std::to_string(d.r->rejoin_cost) + " |";
      } else {
        label_flags += " Dead end rejoin cost=" + std::to_string(d.r->p->cost) +
                       " rjcost=" + std::to_string(d.r->rejoin_cost) + " |";
      }
    }

    if (!d.assembly_ids.empty() || !d.head_assembly_ids.empty()) {
      label_flags += " id=";
      for (auto id : d.head_assembly_ids) {
        label_flags += std::to_string(id) + "(asm) ";
      }
      for (auto id : d.assembly_ids) {
        label_flags += std::to_string(id) + " ";
      }
    }

    if (d.np.new_path->min_overlap != std::numeric_limits<unsigned>::max()) {
      label_flags += " ol:" + std::to_string(d.np.new_path->min_overlap);
    }

    label_flags += " cost:" + std::to_string(d.np.new_path->cost);

    for (const auto& tag : d.tags) {
      label_flags += " " + tag;
    }

#define NP_INFO(FIELD, SHORT)                                                     \
  do {                                                                            \
    if (d.np.FIELD) {                                                             \
      label_flags += " " + std::string(SHORT) + ":" + std::to_string(d.np.FIELD); \
    }                                                                             \
  } while (0)
    NP_INFO(pushed_since_read, "psr");
    NP_INFO(pushed_since_pair, "psp");
    NP_INFO(ambiguous_bases, "a");
    NP_INFO(branch_count_since_pair, "bc");
    NP_INFO(max_between_pairs, "maxbp");
    NP_INFO(pairs_used, "pu");
    NP_INFO(num_reads, "nr");
#undef NP_INFO
    os << label_flags;

    if (!d.costs.empty()) {
      if (!label_flags.empty()) {
        os << " | ";
      }
      os << " | costs:";
      for (const auto& c : d.costs) {
        os << " " << c.first << ":" << c.second;
      }
    }

    os << "}}\"];\n";
  }
  os << "}\n";
}

}  // namespace variants

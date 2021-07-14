#include "modules/variants/filter_dup_align.h"

#include "absl/container/flat_hash_map.h"

namespace variants {

namespace {

constexpr bool k_dbg = false;

}  // namespace

filter_dup_align::filter_dup_align(const sort_func_t& sort_func, pipeline_step_t output)
    : apply_block_step(std::move(output)), m_sort_func(sort_func) {}
filter_dup_align::~filter_dup_align() { flush(); }

struct filter_dup_align::trace_state {
  trace_state(assembly* a_arg) : a(a_arg), filtered_coverage(a->seq.size()) {}

  assembly* const a;
  size_t priority = 0;

  read_id_set read_ids;
  read_coverage_t filtered_coverage;

  std::vector<trace_state*> left_edges;
  std::vector<trace_state*> right_edges;
};

struct filter_dup_align::search_entry {
  aoffset_t cost = 0;
  size_t search_step = 0;

  std::vector<trace_state*> left_path;
  aoffset_t left_offset = 0;
  bool left_end = false;

  std::vector<trace_state*> right_path;
  aoffset_t right_offset = 0;
  bool right_end = false;

  read_coverage_t reads;

  search_entry without_reads() const {
    search_entry result;
    result.cost = cost;
    result.left_path = left_path;
    result.left_offset = left_offset;
    result.left_end = left_end;
    result.right_path = right_path;
    result.right_offset = right_offset;
    result.right_end = right_end;
    return result;
  }
  search_entry() = default;

  // Don't accidentally copy this structure.
  search_entry(const search_entry&) = delete;
  search_entry& operator=(const search_entry&) = delete;

  // But moving is OK.
  search_entry(search_entry&&) = default;
  search_entry& operator=(search_entry&&) = default;

  // Returns true if rhs is preferred in the search.
  bool operator<(const search_entry& rhs) const {
    if (cost != rhs.cost) {
      // Prefer smaller cost.
      return rhs.cost < cost;
    }

    if (left_end != rhs.left_end) {
      // Process ends first.
      return rhs.left_end;
    }

    if (right_end != rhs.right_end) {
      // Process ends first.
      return rhs.right_end;
    }

    // Otherwise, be deterministic by preferring search entries that were created first.
    return rhs.search_step < search_step;
  }
};

class filter_dup_align::block_tracer {
 public:
  block_tracer(aoffset_t left_offset, aoffset_t right_offset, std::vector<assembly_ptr> block)
      : m_block(std::move(block)),
        m_block_left_offset(left_offset),
        m_block_right_offset(right_offset) {}
  void init_state();
  void init_edges();
  void calc_priorities(sort_func_t sort_func);
  void trace_all();
  void flush();

 private:
  trace_state* assembly_to_state(assembly* a);
  void trace(trace_state* st);
  void save_filtered(trace_state* st, const search_entry& entry);
  void add_search_entry(trace_state* st, search_entry entry);
  void search_more(trace_state* st, search_entry entry);
  std::vector<assembly_ptr> m_block;
  std::vector<trace_state*> m_states;
  absl::flat_hash_map<assembly*, std::unique_ptr<trace_state>> m_asm_to_state;

  big_read_id_set m_done_reads;
  const aoffset_t m_block_left_offset;
  const aoffset_t m_block_right_offset;
  absl::btree_multiset<search_entry> m_search_entries;
  size_t m_search_steps = 0;
};

filter_dup_align::trace_state* filter_dup_align::block_tracer::assembly_to_state(assembly* a) {
  auto it = m_asm_to_state.find(a);
  CHECK(it != m_asm_to_state.end());
  return &*it->second;
}

void filter_dup_align::block_tracer::init_state() {
  m_asm_to_state.reserve(m_block.size());

  for (const auto& a : m_block) {
    auto new_state = make_unique<trace_state>(a.get());
    bool did_insert = m_asm_to_state.insert(std::make_pair(a.get(), std::move(new_state))).second;
    CHECK(did_insert) << "Duplicate assembly " << *a << "?";
  }
}

void filter_dup_align::block_tracer::init_edges() {
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

void filter_dup_align::block_tracer::calc_priorities(sort_func_t sort_func) {
  CHECK(m_states.empty());
  size_t old_size = m_block.size();
  m_block = sort_func(std::move(m_block));
  CHECK_EQ(old_size, m_block.size()) << "Sorting should not change assemblies";

  size_t prio = 1;
  for (const auto& a : m_block) {
    auto* st = assembly_to_state(a.get());
    CHECK_EQ(st->priority, 0) << "Duplicate assembly seen when calculating priorities?" << *a;
    st->priority = prio;
    ++prio;
    m_states.push_back(st);
  }
}

void filter_dup_align::block_tracer::trace_all() {
  for (auto* st : m_states) {
    trace(st);
  }
}

void filter_dup_align::block_tracer::flush() {
  for (auto* st : m_states) {
    *st->a->read_coverage = std::move(st->filtered_coverage);
  }
}

void filter_dup_align::block_tracer::trace(trace_state* st) {
  CHECK(m_search_entries.empty());

  big_read_id_set left_to_place;
  {
    search_entry head;
    read_id_set st_reads = st->a->read_coverage->all_read_ids();
    head.reads = *st->a->read_coverage - (st_reads & m_done_reads);
    if (head.reads.empty()) {
      return;
    }
    left_to_place |= head.reads.all_read_ids();
    if (k_dbg) {
      std::cerr << "\nPlacing reads " << head.reads << " in assembly " << *st->a << "\n";
    }
    head.right_offset = st->a->seq.size();
    add_search_entry(st, std::move(head));
  }

  while (!m_search_entries.empty() && !left_to_place.empty()) {
    auto it = m_search_entries.begin();
    search_entry cur = std::move(m_search_entries.extract(it).value());

    read_id_set cur_read_ids = cur.reads.all_read_ids() & left_to_place;
    if (cur_read_ids.empty()) {
      continue;
    }

    cur.reads &= cur_read_ids;
    CHECK(!cur.reads.empty());

    if (k_dbg) {
      std::cerr << "  Placing reads: " << cur.reads << " with " << cur.left_path.size()
                << " left and " << cur.right_path.size() << " right; left_end=" << cur.left_end
                << " right_end=" << cur.right_end << ", search span=[" << cur.left_offset << ", "
                << cur.right_offset << ")\n";
    }

    if (cur.left_end && cur.right_end) {
      if (k_dbg) {
        std::cerr << "   Placing " << cur.reads << "; no more search\n";
      }
      left_to_place -= cur_read_ids;
      save_filtered(st, cur);
      continue;
    }

    search_more(st, std::move(cur));
  }

  m_done_reads |= st->a->read_coverage->all_read_ids();
  m_search_entries.clear();
}

void filter_dup_align::block_tracer::save_filtered(trace_state* st, const search_entry& entry) {
  st->filtered_coverage |= entry.reads;
  aoffset_t left_offset = 0;
  for (const auto& left : entry.left_path) {
    read_coverage_t cov = entry.reads.get_and_adjust_reads_spanning_offset(left_offset);
    aoffset_t seq_size = left->a->seq.size();
    cov.adjust_in_place(seq_size);
    if (k_dbg) {
      std::cerr << "    Placing to the left in assembly id=" << left->a->assembly_id << " ["
                << left->a->left_offset << ", " << left->a->right_offset << "): " << cov << "\n";
    }
    left->filtered_coverage |= cov;
    left_offset -= seq_size;
  }
  CHECK_EQ(left_offset, entry.left_offset);

  aoffset_t right_offset = st->a->seq.size();
  for (const auto& right : entry.right_path) {
    read_coverage_t cov = entry.reads.get_and_adjust_reads_spanning_offset(right_offset);
    aoffset_t seq_size = right->a->seq.size();
    if (k_dbg) {
      std::cerr << "    Placing to the right in assembly id=" << right->a->assembly_id << " ["
                << right->a->left_offset << ", " << right->a->right_offset << "): " << cov << "\n";
    }
    right->filtered_coverage |= cov;
    right_offset += seq_size;
  }
  CHECK_EQ(right_offset, entry.right_offset);
}

void filter_dup_align::block_tracer::add_search_entry(trace_state* st, search_entry entry) {
  CHECK(!entry.reads.empty());
  entry.search_step = m_search_steps++;

  // Don't trace off the end of the block.
  trace_state* right_st = entry.right_path.empty() ? st : entry.right_path.back();
  CHECK_LE(right_st->a->right_offset, m_block_right_offset);
  if (right_st->a->right_offset == m_block_right_offset) {
    entry.right_end = true;
  }

  trace_state* left_st = entry.left_path.empty() ? st : entry.left_path.back();
  CHECK_GE(left_st->a->left_offset, m_block_left_offset);
  if (right_st->a->left_offset == m_block_left_offset) {
    entry.left_end = true;
  }

  m_search_entries.emplace(std::move(entry));
}

void filter_dup_align::block_tracer::search_more(trace_state* st, search_entry entry) {
  if (entry.reads.empty()) {
    return;
  }
  if (!entry.right_end) {
    read_coverage_t past_right = entry.reads.get_reads_spanning_offset(entry.right_offset);
    read_coverage_t right_ends = entry.reads - past_right;
    if (k_dbg) {
      std::cerr << "Searching right past=" << past_right << " ends=" << right_ends << "\n";
    }
    if (!right_ends.empty()) {
      search_entry right_entry = entry.without_reads();
      right_entry.reads = std::move(right_ends);
      right_entry.right_end = true;
      add_search_entry(st, std::move(right_entry));
    }
    if (!past_right.empty()) {
      auto* right_st = entry.right_path.empty() ? st : entry.right_path.back();
      for (trace_state* right_edge : right_st->right_edges) {
        read_coverage_t matching_reads = past_right.intersection_with_adjusted(
            *right_edge->a->read_coverage, entry.right_offset);
        if (matching_reads.empty()) {
          continue;
        }
        search_entry right_entry = entry.without_reads();
        right_entry.reads = std::move(matching_reads);

        right_entry.cost += right_edge->priority;
        right_entry.right_path.push_back(right_edge);
        right_entry.right_offset += right_edge->a->seq.size();
        add_search_entry(st, std::move(right_entry));
      }
    }
  }

  if (!entry.left_end) {
    read_coverage_t past_left = entry.reads.get_reads_spanning_offset(entry.left_offset);
    read_coverage_t left_ends = entry.reads - past_left;
    if (k_dbg) {
      std::cerr << "Searching left past=" << past_left << " ends=" << left_ends << "\n";
    }
    if (!left_ends.empty()) {
      search_entry left_entry = entry.without_reads();
      left_entry.reads = std::move(left_ends);
      left_entry.left_end = true;
      add_search_entry(st, std::move(left_entry));
    }
    if (!past_left.empty()) {
      auto* left_st = entry.left_path.empty() ? st : entry.left_path.back();
      for (trace_state* left_edge : left_st->left_edges) {
        read_coverage_t matching_reads = past_left.intersection_with_adjusted(
            *left_edge->a->read_coverage, entry.left_offset - aoffset_t(left_edge->a->seq.size()));
        if (matching_reads.empty()) {
          continue;
        }
        search_entry left_entry = entry.without_reads();
        left_entry.reads = std::move(matching_reads);
        left_entry.cost += left_edge->priority;
        left_entry.left_path.push_back(left_edge);
        left_entry.left_offset -= aoffset_t(left_edge->a->seq.size());
        add_search_entry(st, std::move(left_entry));
      }
    }
  }
}

void filter_dup_align::on_block(aoffset_t left_offset, aoffset_t right_offset,
                                const std::vector<assembly_ptr>& block) {
  if (k_dbg) {
    std::cerr
        << "********************************************************************************\n";
    std::cerr << "Got block [" << left_offset << ", " << right_offset << ") with " << block.size()
              << " assemblies\n";
  }
  if (block.size() <= 1) {
    // No need to resolve anything.
    return;
  }

  block_tracer tracer(left_offset, right_offset, block);
  tracer.init_state();
  tracer.init_edges();
  tracer.calc_priorities(m_sort_func);
  tracer.trace_all();
  tracer.flush();
}

}  // namespace variants

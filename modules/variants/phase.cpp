#include "modules/variants/phase.h"
#include "modules/io/ref_count.h"
#include "modules/io/stats.h"

namespace variants {

namespace {

constexpr bool k_dbg = false;
constexpr bool k_show_stats = false;
constexpr int k_show_stats_interval = 30;

static time_t g_next_show_stats = 0;
static time_t g_next_split_show_stats = 0;

}  // namespace

const char join_phases::k_join_phases_name[] = "JOIN_PHASES";

bool join_phases::g_check_invariants = false;

join_phases::join_phases(size_t max_phase_len, size_t max_phase_asm_len, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)),
      m_max_phase_len(max_phase_len),
      m_max_phase_asm_len(max_phase_asm_len) {}

join_phases::~join_phases() { flush(); }

void join_phases::show_stats() {
  if (!k_show_stats) {
    return;
  }

  time_t now = time(0);
  if (g_next_show_stats > now) {
    return;
  }
  g_next_show_stats += k_show_stats_interval;
  if (g_next_show_stats < now) {
    g_next_show_stats = now + k_show_stats_interval;
    if (g_next_split_show_stats < now) {
      g_next_split_show_stats = now;
    }
  }

  absl::btree_set<active_t*> distinct_actives;
  absl::btree_set<const assembly*> distinct_subs;
  absl::btree_set<const assembly*> distinct_refs;

  for (const auto& act : m_active) {
    distinct_actives.insert(act.second.get());
  }
  simple_stats<double> dist_to_left;
  simple_stats<double> dist_to_right;
  simple_stats<double> dist_to_var_right;
  simple_stats<double> num_joined;
  simple_stats<double> num_ref;
  size_t tot_joined = 0;
  size_t tot_ref = 0;
  for (active_t* act : distinct_actives) {
    dist_to_left.add_sample(act->joined_a->left_offset - m_cur_offset);
    dist_to_right.add_sample(act->right_offset - m_cur_offset);
    dist_to_var_right.add_sample(act->var_right_offset - m_cur_offset);
    num_joined.add_sample(act->joined_a->sub_assemblies.size());
    tot_joined += act->joined_a->sub_assemblies.size();
    num_ref.add_sample(act->reference_after.size());
    tot_ref += act->reference_after.size();

    for (const auto& sub : act->joined_a->sub_assemblies) {
      distinct_subs.insert(sub->get());
    }
    for (const auto& ref : act->reference_after) {
      distinct_refs.insert(ref->get());
    }
  }

  auto statstr = [](simple_stats<double>& stats) -> std::string {
    std::stringstream out;
    stats.analyze();
    if (stats.samples.empty()) {
      out << "(none)";
    } else {
      out << stats.min << "/" << stats.avg << "/" << stats.max;
    }
    return out.str();
  };

  std::cerr << "join_phases@" << m_cur_offset << ": " << m_tot_seen << " asms (" << tot_joined
            << " active, " << distinct_subs.size() << " distinct), " << m_tot_seen_phases
            << " phases (" << m_seen_phases.size() << " distinct, " << m_active.size()
            << " active, " << distinct_actives.size() << " distinct active) ref afters=" << tot_ref
            << " (" << distinct_refs.size() << " distinct)"
            << " left=" << statstr(dist_to_left) << " right=" << statstr(dist_to_right)
            << " var_right=" << statstr(dist_to_var_right) << "\n";
}

void join_phases::check_invariants() {
  if (!g_check_invariants) {
    return;
  }

  absl::btree_map<active_t*, phase_set> phase_refs;
  for (const auto& act : m_active) {
    CHECK(phase_refs[act.second.get()].insert(act.first).second);
  }

  for (const auto& phase_ref : phase_refs) {
    active_t* act = phase_ref.first;
    CHECK(phase_ref.second == act->joined_a->phase_ids)
        << "phase_ref: " << phase_ref.second << " active: " << act->joined_a->phase_ids;

    aoffset_t right_offset = act->joined_a->left_offset;
    dna_sequence seq;
    for (const auto& sub : act->joined_a->sub_assemblies) {
      CHECK_EQ((*sub)->left_offset, right_offset);
      right_offset = (*sub)->right_offset;
      seq += (*sub)->seq;
    }
    CHECK_EQ(act->joined_a->right_offset, right_offset);
    CHECK_EQ(seq, act->joined_a->seq);
    CHECK_EQ(right_offset, act->var_right_offset);

    for (const auto& ref_asm : act->reference_after) {
      CHECK_EQ((*ref_asm)->left_offset, right_offset);
      right_offset = (*ref_asm)->right_offset;
    }
    CHECK_EQ(right_offset, act->right_offset);
  }
}

void join_phases::add_ref_asm(std::shared_ptr<assembly_ptr> shared_a) {
  if (k_dbg) {
    std::cout << "join_phases::add_ref_asm: " << **shared_a << "\n";
  }
  CHECK((*shared_a)->matches_reference);
  auto& a = *shared_a;

  for (auto it = m_active.begin(); it != m_active.end();) {
    auto& act = it->second;

    if (act->right_offset > a->left_offset) {
      if (k_dbg) {
        std::cout << "Phase " << it->first << " at " << act->right_offset
                  << "; doesn't need reference at " << a->left_offset << "\n";
      }
      ++it;
      continue;
    }

    if (act->right_offset < a->left_offset) {
      if (k_dbg) {
        std::cout << "Phase " << it->first << " at " << act->right_offset << "; not caught up to "
                  << a->left_offset << ". discarding.\n";
      }
      output_active(std::move(act));
      it = m_active.erase(it);
      continue;
    }

    if (k_dbg) {
      std::cout << "Phase " << it->first << " at " << act->right_offset << " adding reference from "
                << a->left_offset << " to " << a->right_offset << "\n";
    }
    CHECK_EQ((*shared_a)->left_offset, act->right_offset);
    act->reference_after.push_back(shared_a);
    act->right_offset = (*shared_a)->right_offset;
    if (k_dbg) {
      std::cout << "new active right offset: " << act->right_offset << "\n";
    }
    ++it;
  }
  output_if_last_ref(std::move(shared_a));

  check_invariants();
}

void join_phases::output_active(active_ptr act_shared) {
  auto act = act_shared.release();
  if (!act) {
    // Not the last reference; wait to output it.
    return;
  }
  if (k_dbg) {
    std::cout << "Outputting active " << *act->joined_a << "\n";
  }
  untrack_left_offset(act->joined_a->left_offset);
  sort_and_output(std::move(act->joined_a));
  for (auto& ref : act->reference_after) {
    output_if_last_ref(std::move(ref));
  }
}

void join_phases::output_if_last_ref(std::shared_ptr<assembly_ptr> shared_a) {
  if (shared_a.use_count() == 1) {
    assembly_ptr a = std::move(*shared_a);
    if (k_dbg) {
      std::cout << "Outputting since it's the last reference: " << *a << "\n";
    }
    sort_and_output(std::move(a));
  } else {
    if (k_dbg) {
      std::cout << shared_a.use_count() << " refs left on " << **shared_a << "\n";
    }
  }
}

void join_phases::advance_to(aoffset_t target) {
  if (k_dbg) {
    std::cout << "Advance target: " << target << "\n";
  }
  while (m_cur_offset < target) {
    advance_towards(target);
    flush_sorted_to(m_cur_offset);
    show_stats();
  }
  CHECK_EQ(m_cur_offset, target);
}

void join_phases::advance_towards(aoffset_t target) {
  CHECK_GT(target, m_cur_offset);
  for (auto& ref_asm : m_cur_ref) {
    aoffset_t left_offset = ref_asm->left_offset;

    auto shared_a = std::make_shared<assembly_ptr>(std::move(ref_asm));
    add_ref_asm(shared_a);

    auto new_act = new_active((*shared_a)->left_offset, {});
    CHECK_EQ((*shared_a)->left_offset, new_act->right_offset);
    new_act->right_offset = (*shared_a)->right_offset;
    new_act->joined_a->matches_reference = true;
    add_to_active(new_act.get(), std::move(shared_a));

    CHECK_EQ(1, new_act.use_count()) << "Reference phase should not be added to any phase group";
    output_active(std::move(new_act));
    untrack_left_offset(left_offset);
  }
  m_cur_ref.clear();

  for (const auto& act_elem : m_active) {
    const auto& act = act_elem.second;
    if (act->right_offset > m_cur_offset) {
      target = std::min(target, act->right_offset);
    }
  }

  if (!m_abort_at.empty()) {
    target = std::min(target, m_abort_at.begin()->first);
  }

  CHECK_GT(target, m_cur_offset);
  if (k_dbg) {
    std::cout << "Advancing " << m_cur_offset << " to " << target << "\n";
  }
  m_cur_offset = target;

  for (auto it = m_active.begin(); it != m_active.end();) {
    auto& act = it->second;

    if (act->right_offset < m_cur_offset ||
        act->var_right_offset + aoffset_t(m_max_phase_len) < m_cur_offset) {
      if (k_dbg) {
        std::cout << act->right_offset << " (var=" << act->var_right_offset
                  << ") expired phase: " << *act->joined_a << " ids=" << act->joined_a->phase_ids
                  << "\n";
      }
      output_active(std::move(act));
      it = m_active.erase(it);

      continue;
    }

    ++it;
  }

  if (!m_abort_at.empty()) {
    auto first_abort = m_abort_at.begin();
    CHECK_GE(first_abort->first, m_cur_offset);
    if (first_abort->first == m_cur_offset) {
      abort_phases(first_abort->second);
      m_abort_at.erase(first_abort);
    }
  }
  check_invariants();
}

void join_phases::on_assembly(assembly_ptr a) {
  if (k_dbg) {
    std::cout << "join_phases got assembly " << *a << " with phases " << a->phase_ids << "\n";
  }
  if (k_show_stats) {
    ++m_tot_seen;
    m_tot_seen_phases += a->phase_ids.size();
    for (const std::string& phase_id : a->phase_ids) {
      m_seen_phases.insert(phase_id);
    }
  }
  advance_to(a->left_offset);

  if (a->phase_ids.empty() && !a->matches_reference) {
    if (k_dbg) {
      std::cout << "Unphased assembly passing through: " << *a << "\n";
    }
    sort_and_output(std::move(a));
    return;
  }

  if (a->matches_reference) {
    if (k_dbg) {
      std::cout << "Saving ref asm: " << *a << "\n";
    }
    // Process reference last, after all variants here.
    track_left_offset(a->left_offset);
    m_cur_ref.emplace_back(std::move(a));
    return;
  }

  add_var_asm(std::make_shared<assembly_ptr>(std::move(a)));
}

void join_phases::add_var_asm(std::shared_ptr<assembly_ptr> shared_a) {
  if (k_dbg) {
    std::cout << "Adding var asm: " << **shared_a << "\n";
  }

  auto& a = *shared_a;
  CHECK(!a->phase_ids.empty());

  bool force_abort = false;
  if ((a->right_offset - a->left_offset) > aoffset_t(m_max_phase_asm_len)) {
    if (k_dbg) {
      std::cout << "Forcing abort due to large variant seq\n";
    }
    force_abort = true;
  } else if (aoffset_t(a->seq.size()) > m_max_phase_asm_len) {
    if (k_dbg) {
      std::cout << "Forcing abort due to large reference seq\n";
    }
    force_abort = true;
  }

  absl::btree_map<active_t*, phase_set /* phase ids */> found_phases;
  absl::btree_map<active_t*, std::pair<active_ptr, phase_set /* phase ids */>> discard_phases;
  phase_set new_phases;
  phase_set abort_phase_ids;
  if (force_abort) {
    abort_phase_ids = a->phase_ids;
  } else {
    for (const auto& phase_id : a->phase_ids) {
      auto it = m_active.find(phase_id);
      if (it == m_active.end()) {
        new_phases.insert(phase_id);
        continue;
      }
      auto& act = it->second;
      if (act->right_offset != a->left_offset) {
        std::stringstream msg;
        msg << "Phase conflict with phase id '" << phase_id << "' between " << *act->joined_a
            << " and " << *a << "; consider running resolve_phase_conflicts";
        throw(io_exception(msg.str()));
      }

      found_phases[act.get()].insert(phase_id);
    }
  }

  if (!abort_phase_ids.empty()) {
    abort_phases(abort_phase_ids);
    // Abort these phase ids again at the other end of this assembly.
    if (a->right_offset > m_cur_offset) {
      m_abort_at[a->right_offset] += abort_phase_ids;
    }

    new_phases += abort_phase_ids;
  }

  for (const auto& found : found_phases) {
    auto* act = found.first;
    const auto& found_phase_ids = found.second;

    if (k_dbg) {
      std::cout << "Processing phase subset: " << found_phase_ids << "\n";
    }

    if (act->joined_a->phase_ids != found_phase_ids) {
      // Split necessary.
      split_active(act, found_phase_ids);
      CHECK(act->joined_a->phase_ids == found_phase_ids)
          << "active=" << act->joined_a->phase_ids << " found=" << found_phase_ids;
    }

    CHECK_EQ(act->right_offset, a->left_offset);
    save_ref_asms(act);
    add_to_active(act, shared_a);
    CHECK_EQ((*shared_a)->left_offset, act->right_offset);
    act->right_offset = (*shared_a)->right_offset;
    check_invariants();
  }

  if (new_phases.empty()) {
    CHECK(!found_phases.empty());
  } else {
    auto new_act = new_active((*shared_a)->left_offset, new_phases);
    CHECK_EQ((*shared_a)->left_offset, new_act->right_offset);
    new_act->right_offset = (*shared_a)->right_offset;
    add_to_active(new_act.get(), std::move(shared_a));
    CHECK(!new_act.release()) << "Unable to add new active to any phases";
  }

  if (!abort_phase_ids.empty()) {
    abort_phases(abort_phase_ids);
  }
  check_invariants();
}

void join_phases::abort_phases(const phase_set& abort_ids) {
  if (k_dbg) {
    std::cout << "Aborting phases: " << abort_ids << "\n";
  }

  CHECK(!abort_ids.empty());
  absl::btree_map<active_t*, std::pair<active_ptr, phase_set /* phase ids */>> found_phases;
  for (const auto& abort_id : abort_ids) {
    auto it = m_active.find(abort_id);
    if (it == m_active.end()) {
      // No need to abort this phase id.
      continue;
    }
    auto& act = it->second;

    auto& found = found_phases[act.get()];
    auto& found_act = found.first;
    auto& found_phase_ids = found.second;

    if (found_act) {
      CHECK_EQ(found_act.get(), act.get());
      CHECK(!found_act.release());
    }
    found_act = std::move(act);
    CHECK(found_phase_ids.insert(abort_id).second);

    m_active.erase(it);
  }

  for (auto& found : found_phases) {
    auto& act = found.second.first;
    CHECK_EQ(found.first, act.get());
    auto& phase_ids = found.second.second;

    if (act->joined_a->phase_ids != phase_ids) {
      split_active(act.get(), phase_ids);
    }
    CHECK_EQ(act->joined_a->phase_ids, phase_ids);
    output_active(std::move(act));
  }
}

void join_phases::split_active(active_t* act, const phase_set& keep_phases) {
  CHECK_LT(keep_phases.size(), act->joined_a->phase_ids.size());

  if (g_check_invariants) {
    phase_set unexpected_phases = keep_phases - act->joined_a->phase_ids;
    CHECK(unexpected_phases.empty())
        << "orig=" << act->joined_a->phase_ids << "keep=" << keep_phases
        << "unexpected=" << unexpected_phases;
  }

  phase_set split_phases = act->joined_a->phase_ids - keep_phases;
  CHECK(!split_phases.empty()) << "orig=" << act->joined_a->phase_ids << "keep=" << keep_phases
                               << "split=" << split_phases;

  if (k_dbg) {
    std::cout << "Splitting active " << *act->joined_a << ", keeping " << keep_phases
              << ", splitting " << split_phases << "\n";
  }

  act->joined_a->phase_ids = keep_phases;
  active_ptr new_act(make_unique<active_t>());
  new_act->joined_a = make_unique<assembly>(*act->joined_a);
  new_act->joined_a->assembly_id = allocate_assembly_id();
  new_act->joined_a->phase_ids = std::move(split_phases);
  new_act->reference_after = act->reference_after;
  new_act->right_offset = act->right_offset;
  new_act->var_right_offset = act->var_right_offset;
  track_left_offset(new_act->joined_a->left_offset);

  for (const auto& phase_id : new_act->joined_a->phase_ids) {
    auto it = m_active.find(phase_id);
    CHECK(it != m_active.end());
    CHECK_EQ(it->second.get(), act);
    CHECK(!it->second.release()) << "Split removed last old active?";
    it->second = new_act.clone();
  }
  CHECK(!new_act.release()) << "Unable to assign new active to any phase ids";
  check_invariants();
}

void join_phases::save_ref_asms(active_t* act) {
  if (k_dbg) {
    std::cout << "Saving " << act->reference_after.size() << " ref asms from " << *act->joined_a
              << " init right = " << act->right_offset << "\n";
  }

  for (auto& ref_asm : act->reference_after) {
    add_to_active(act, std::move(ref_asm));
  }
  act->reference_after.clear();
}

void join_phases::add_to_active(active_t* act, std::shared_ptr<assembly_ptr> shared_a) {
  const auto& a = *shared_a;
  act->joined_a->seq += a->seq;
  CHECK_EQ(act->var_right_offset, a->left_offset);
  CHECK_EQ(act->joined_a->right_offset, a->left_offset);
  act->joined_a->right_offset = a->right_offset;
  act->joined_a->sub_assemblies.push_back(std::move(shared_a));
  act->var_right_offset = a->right_offset;
}

join_phases::active_ptr join_phases::new_active(aoffset_t left_offset,
                                                const phase_set& new_phases) {
  if (k_dbg) {
    std::cout << "New active for " << new_phases << " at " << left_offset << "\n";
  }

  active_ptr new_act(make_unique<active_t>());
  new_act->joined_a = make_unique<assembly>();
  new_act->right_offset = left_offset;
  new_act->var_right_offset = left_offset;
  auto& new_a = new_act->joined_a;
  new_a->assembly_id = allocate_assembly_id();
  new_a->left_offset = left_offset;
  new_a->right_offset = left_offset;
  new_a->tags.insert(k_join_phases_name);
  new_a->phase_ids = new_phases;
  track_left_offset(new_a->left_offset);

  for (const auto& phase_id : new_phases) {
    auto& act = m_active[phase_id];
    CHECK(!act);
    act = new_act.clone();
  }
  return new_act;
}

void join_phases::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
  CHECK(m_cur_ref.empty());
}

void propagate_subassembly_coverage(assembly* a) {
  if (a->sub_assemblies.empty()) {
    return;
  }

  aoffset_t seq_offset = 0;

  for (auto& subasm : a->sub_assemblies) {
    const auto& suba = *subasm;

    size_t subseq_size = (*subasm)->seq.size();

    auto propagate_cov = [seq_offset, subseq_size](boost::optional<read_coverage_t>& cov,
                                                   boost::optional<read_coverage_t>& subcov) {
      if (!cov) {
        return;
      }

      if (subcov) {
        // subcov->add_subcoverage(*cov, seq_offset, subseq_size);
        read_coverage_t new_subcov = cov->subcoverage(seq_offset, subseq_size);
        *subcov = subcov->union_with(new_subcov);
      } else {
        read_coverage_t new_subcov = cov->subcoverage(seq_offset, subseq_size);
        subcov.emplace(std::move(new_subcov));
      }
    };
    propagate_cov(a->read_coverage, suba->read_coverage);
    propagate_cov(a->pair_read_coverage, suba->pair_read_coverage);
    seq_offset += subseq_size;
  }
  CHECK_EQ(seq_offset, a->seq.size());
}

split_phases::split_phases(pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)) {}

split_phases::~split_phases() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
}

void split_phases::on_assembly(assembly_ptr a) {
  advance_to(a->left_offset);
  if (k_show_stats) {
    time_t now = time(0);
    if (g_next_split_show_stats < now) {
      g_next_split_show_stats += k_show_stats_interval;
      if (g_next_split_show_stats < now) {
        g_next_split_show_stats = now + k_show_stats_interval;
      }
      std::cerr << "split_phases@" << a->left_offset << " active=" << m_active.size();
      if (!m_active.empty()) {
        std::cerr << " left=" << m_active.begin()->first << " right=" << m_active.rbegin()->first;
      }
      std::cerr << "\n";
    }
  }
  if (a->sub_assemblies.empty()) {
    sort_and_output(std::move(a));
    return;
  }

  for (auto& suba : a->sub_assemblies) {
    aoffset_t left_offset = (*suba)->left_offset;
    if (m_active[left_offset].insert(std::move(suba)).second) {
      track_left_offset(left_offset);
    }
  }
}

void split_phases::advance_to(aoffset_t offset) {
  while (!m_active.empty() && m_active.begin()->first < offset) {
    auto& asms = m_active.begin()->second;
    for (const auto& a : asms) {
      untrack_left_offset((*a)->left_offset);
      sort_and_output(std::move(*a));
    }
    m_active.erase(m_active.begin());
  }
  flush_sorted_to(offset);
}

resolve_phase_conflicts::resolve_phase_conflicts(const resolve_conflict_func_t& resolve_conflict,
                                                 pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_resolve_conflict(resolve_conflict) {}

resolve_phase_conflicts::~resolve_phase_conflicts() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
}

void resolve_phase_conflicts::on_assembly(assembly_ptr a) {
  advance_to(a->left_offset);

  bool need_conflict_check = true;
  while (need_conflict_check) {
    need_conflict_check = false;
    for (auto& act : m_active) {
      if (check_and_resolve_conflicts(a, act.second)) {
        need_conflict_check = true;
        break;
      }
    }
  }

  track_left_offset(a->left_offset);
  aoffset_t right_offset = a->right_offset;
  m_active.emplace(right_offset, std::move(a));
}

void resolve_phase_conflicts::advance_to(aoffset_t pos) {
  while (!m_active.empty() && m_active.begin()->first < pos) {
    auto it = m_active.begin();
    assembly_ptr a = std::move(it->second);
    m_active.erase(it);
    untrack_left_offset(a->left_offset);
    sort_and_output(std::move(a));
  }
  flush_sorted_to(pos);
}

bool resolve_phase_conflicts::check_and_resolve_conflicts(const assembly_ptr& a,
                                                          const assembly_ptr& b) {
  phase_set in_common = a->phase_ids & b->phase_ids;
  if (in_common.empty()) {
    return false;
  }

  m_resolve_conflict(a, b, in_common);

  phase_set new_in_common = a->phase_ids & b->phase_ids;
  if (!new_in_common.empty()) {
    std::stringstream msg;
    msg << "Phase conflict resolution failed between " << *a << " and " << *b
        << "; origin conflict phases=" << in_common
        << ", after resolution conflicts=" << new_in_common;
    throw(io_exception(msg.str()));
  }

  return true;
}

}  // namespace variants

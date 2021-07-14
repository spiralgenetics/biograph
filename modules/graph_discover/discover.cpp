#include "modules/graph_discover/discover.h"

#include <boost/range/adaptor/reversed.hpp>
#include <chrono>
#include <ctime>
#include <random>

#include "modules/bio_base/readmap.h"
#include "modules/io/parallel.h"

namespace variants {

constexpr int k_dbg = 0;

std::ostream& operator<<(std::ostream& os, const graph_discover::active_assembly& act) {
  return os << act.to_string();
}

std::string graph_discover::active_assembly::to_string() const {
  std::stringstream os;
  os << "Graph discover active asm=";
  if (a) {
    os << *a << "\n";
  } else {
    os << "(no assembly)\n";
  }

  return os.str();
}

graph_discover::graph_discover(const assemble_options& options, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_options(options) {
  CHECK(m_options.seqset);
  CHECK(m_options.readmap);
}

void graph_discover::flush() {
  advance_trace_to(std::numeric_limits<aoffset_t>::max());
  flush_sorted();
}

graph_discover::~graph_discover() {
  flush();
  CHECK_EQ(m_trace_offset, std::numeric_limits<aoffset_t>::max());
}

void graph_discover::process_readahead(const active_assembly* act) {
  if (k_dbg) {
    std::cerr << "Readahead: " << *act << "\n";
  }
  on_readahead(act);
}

void graph_discover::process_readahead_done(const active_assembly* act) {
  // Remove references to this assembly.

  if (k_dbg) {
    std::cerr << "Readahead done: " << *act << "\n";
  }

  on_readahead_done(act);
}

void graph_discover::process_trace(const active_assembly* act) {
  if (!opts().discover_tags.empty()) {
    if ((opts().discover_tags & act->a->tags).empty()) {
      if (k_dbg) {
        std::cerr << "Not tracing due to no matching in common: " << *act << "\n";
      }
      return;
    }
    if (k_dbg) {
      std::cerr << "Tags in common; executing trace: " << *act << "\n";
    }
  } else {
    if (k_dbg) {
      std::cerr << "No tags configured; Executing trace: " << *act << "\n";
    }
  }

  on_trace(act);
}

void graph_discover::on_assembly(assembly_ptr a) {
  track_left_offset(min(a->left_offset, a->right_offset));
  advance_trace_to(min(a->left_offset, a->right_offset) - opts().read_ahead_distance);

  active_assembly_ptr new_act = make_active_assembly();

  seqset_range_set rc_entry_ends;

  if (k_dbg) {
    std::cerr << "On assembly: " << *a << "\n";
  }
  new_act->a = std::move(a);
  walk_readahead(std::move(new_act));
}

void graph_discover::walk_readahead(active_assembly_ptr act) {
  on_walk(act.get());

  process_readahead(act.get());
  optional_aoffset min_offset = min(act->a->left_offset, act->a->right_offset);
  m_readahead_done.emplace(min_offset, std::move(act));
}

void graph_discover::advance_trace_to(aoffset_t pos) {
  while (m_trace_offset < pos) {
    flush_sorted_to(m_trace_offset);
    advance_trace_towards(pos);
  }
}

void graph_discover::advance_trace_towards(aoffset_t pos) {
  if (!m_trace_pending.empty() && m_trace_pending.begin()->first < pos) {
    pos = m_trace_pending.begin()->first;
  }

  if (k_dbg) {
    std::cerr << "Flushing readahead to " << pos << "\n";
  }
  while (!m_readahead_done.empty() && m_readahead_done.begin()->first < pos) {
    active_assembly_ptr act = std::move(m_readahead_done.begin()->second);
    m_readahead_done.erase(m_readahead_done.begin());

    process_readahead_done(act.get());

    aoffset_t max_offset = max(act->a->left_offset, act->a->right_offset);
    m_trace_pending.emplace(max_offset, std::move(act));
    if (max_offset < pos) {
      pos = max_offset;
    }
  }

  if (k_dbg) {
    std::cerr << "Advancing trace from " << m_trace_offset << " to " << pos << "\n";
  }

  m_trace_offset = pos;
  on_advance_trace(m_trace_offset);

  while (!m_trace_pending.empty() && m_trace_pending.begin()->first == pos) {
    active_assembly_ptr act = std::move(m_trace_pending.begin()->second);
    m_trace_pending.erase(m_trace_pending.begin());

    process_trace(act.get());

    untrack_left_offset(min(act->a->left_offset, act->a->right_offset));
    sort_and_output(std::move(act->a));
  }
}

void graph_discover::extend_assembly(const assembly* orig, assembly* a, dna_slice seq) const {
  a->seq += seq;
}

assembly_ptr graph_discover::discover_extend_right(const active_assembly* act, aoffset_t offset,
                                                   dna_slice seq, const std::string& tag,
                                                   seqset_path new_rc_path) {
  assembly_ptr a = make_unique<assembly>();
  a->assembly_id = allocate_assembly_id();
  a->tags.insert(tag);
  a->seq = act->a->seq.subseq(0, offset);
  a->rc_seqset_entries = new_rc_path;
  a->left_offset = act->a->left_offset;
  CHECK(act->a->left_offset);
  a->right_offset = optional_aoffset::none;
  a->matches_reference = false;

  extend_assembly(act->a.get(), a.get(), seq);

  return a;
}

assembly_ptr graph_discover::discover_anchor(const active_assembly* act, aoffset_t offset,
                                             dna_slice seq, const potential_anchor& anchor,
                                             const std::string& tag, seqset_path new_rc_path) {
  assembly_ptr a = make_unique<assembly>();
  a->assembly_id = allocate_assembly_id();
  a->tags.insert(tag);
  a->left_offset = act->a->left_offset;
  CHECK(act->a->left_offset);
  a->seq = act->a->seq.subseq(0, offset);
  a->rc_seqset_entries = new_rc_path;

  CHECK(anchor.act->a->right_offset);
  a->right_offset = anchor.act->a->right_offset;
  a->matches_reference = false;

  dna_sequence ext_seq;
  ext_seq += seq;
  ext_seq += anchor.act->a->seq.subseq(anchor.offset, anchor.act->a->seq.size() - anchor.offset);

  extend_assembly(act->a.get(), a.get(), ext_seq);

  return a;
}

std::ostream& operator<<(std::ostream& os, const seqset_range_set& rs) {
  if (rs.empty()) {
    return os << "(empty)";
  }
  for (const auto& r : rs) {
    os << r.sequence() << " ";
  }
  return os << "(" << rs.size() << " ranges)";
}

}  // namespace variants

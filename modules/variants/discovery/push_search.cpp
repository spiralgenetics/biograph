#include "modules/variants/discovery/push_search.h"
#include "modules/variants/discovery/branch.h"
#include "modules/variants/discovery/path.h"
#include "modules/variants/discovery/pop_search.h"
#include "modules/variants/discovery/rejoin.h"

namespace variants {
namespace discovery {

constexpr bool k_trace_all = false;
constexpr bool k_trace_selected = true;

class push_tracer {
 public:
  push_tracer(push_search_entry* e, branch* br)
      : m_e(e),
        m_br(br),
        m_st(m_br->get_state()),
        m_view(m_br->push_view()),
        m_rc_view(m_br->pop_view()),
        m_opts(m_br->push_view()->opts()),
        m_path(e->m_path) {}

  search_result search() {
    if (k_trace_selected && m_br->trace_enabled(m_path)) {
      m_trace = true;
    }
    if (m_trace) {
      std::cout << "Starting push trace: " << m_e->describe(m_br) << "\n";
    }
    populate_cur_read();
    save_pop_and_rp();
    if (!m_e->m_did_rejoin) {
      check_rejoin();
      if (!m_br->explore(m_path.range())) {
        if (m_trace) {
          std::cout
              << "Already explored this seqset range; don't do this homologous region again.\n";
        }
        return search_result::STOP_SEARCHING;
      }
      flush_last_rp();
      flush_last_pop();
    }

    int orig_cur_overlap = std::min<int>(m_path.last_overlap(), m_path.cur_overlap());
    unsigned orig_path_overlap = m_e->m_key.path_overlap;
    while (trace_one_base()) {
      CHECK(!m_e->m_did_rejoin);
      if (m_e->m_key.path_overlap != orig_path_overlap ||
          m_path.last_overlap() < orig_cur_overlap) {
        CHECK_LE(m_e->m_key.path_overlap, orig_path_overlap)
            << "Path overlap should only decrease while tracing";
        flush_last_rp();
        flush_last_pop();
        if (m_trace) {
          std::cout << "Path overlap decreased from " << orig_path_overlap << " to "
                    << m_e->m_key.path_overlap << ".  Cur overlap decreased from "
                    << orig_cur_overlap << " to " << m_path.last_overlap()
                    << ".  Search more later.\n";
        }
        // Our overlap decreased; continue searching later.
        return search_result::SEARCH_MORE;
      }

      populate_cur_read();
      check_pair_support();
      save_pop_and_rp();
      check_rejoin();
      if (m_ambig) {
        if (m_trace) {
          std::cout << "Ambiguous; stopping search\n";
        }
        flush_last_rp();
        flush_last_pop();
        break;
      }
      if (m_opts.bidir_pop_all_reads) {
        flush_last_pop();
      }
      if (m_opts.bidir_right_partial_all_reads) {
        flush_last_rp();
      }

      if (m_e->m_did_rejoin) {
        if (m_trace) {
          std::cout << "Pausing push search for now because we rejoined\n";
        }
        flush_last_rp();
        return search_result::SEARCH_MORE;
      }

      if (m_got_pair_support_or_mateless) {
        flush_last_pop();
      }
    }
    flush_last_rp();
    flush_last_pop();
    if (m_trace) {
      std::cout << "Stopping push search.\n";
    }
    return search_result::STOP_SEARCHING;
  }

 private:
  void flush_last_rp() {
    if (m_pending_rp) {
      if (m_trace) {
        std::cout << "Flushing right partial: " << *m_pending_rp << "\n";
      }
      m_view->add_right_partial(*m_pending_rp_r, std::move(*m_pending_rp));
      m_pending_rp.reset();
      m_pending_rp_r.reset();
    }
  }

  void flush_last_pop() {
    if (m_pending_pop) {
      if (m_trace) {
        std::cout << "Flushing pop: " << m_pending_pop->describe(m_br) << "\n";
      }
      m_br->add_search_entry(std::move(m_pending_pop));
      m_pending_pop.reset();
    }
  }

  void populate_cur_read() {
    if (!m_seen.insert(m_path.range()).second) {
      if (m_trace) {
        std::cout << "Already explored; setting ambiguous\n";
      }
      m_ambig = true;
    }

    if (!m_path.longest_read_id()) {
      m_cur_read_r = seqset_range();
      return;
    }

    m_cur_read_r = m_path.range();
    uint32_t read_id = *m_path.longest_read_id();
    uint32_t rc_read_id = m_opts.readmap->get_rev_comp(read_id);
    uint64_t seqset_id = m_opts.readmap->index_to_entry(rc_read_id);
    m_cur_read_rc_r =
        m_opts.seqset->ctx_entry(seqset_id).truncate(m_opts.readmap->get_readlength(rc_read_id));

    if (m_trace) {
      std::cout << "Found read: " << m_cur_read_r.sequence() << "\n";
      std::cout << "Path so far: " << m_path << "\n";
    }
  }

  void save_pop_and_rp() {
    if (!m_cur_read_r.valid()) {
      return;
    }

    if (m_got_pair_support_or_mateless) {
      m_pending_pop =
          make_unique<pop_search_entry>(m_cur_read_rc_r, m_path, m_e->pair_match_count());

      aoffset_t outer_right_offset = m_br->right_push_view_offset() + m_path.anchor_len();
      m_pending_rp.emplace(m_path.seq(), outer_right_offset, m_e->pair_match_count());
      m_pending_rp_r.emplace(m_cur_read_r);
    }
  }

  void check_pair_support() {
    if (!m_cur_read_r.valid()) {
      return;
    }
    auto rc_ri_entries = m_rc_view->range_info().entries_starting_with(m_cur_read_rc_r);
    bool got_pair_support = false;
    for (const auto& r_and_ri : rc_ri_entries) {
      const auto& rc_r = r_and_ri.first;
      const auto& rc_ri = r_and_ri.second;
      CHECK_GE(rc_r.begin(), m_cur_read_rc_r.begin());
      CHECK_LE(rc_r.end(), m_cur_read_rc_r.end());
      interval_set_t relevant_rc_supported_offsets =
          rc_ri.pair_supported_offsets &
          interval_t(m_rc_view->reverse_offset(m_br->right_push_view_offset()),
                     m_rc_view->reverse_offset(m_br->push_view_farthest_left_offset()));
      if (!relevant_rc_supported_offsets.empty()) {
        if (m_trace) {
          std::cout << "Push trace found read with pair support\n";
        }
        got_pair_support = true;
      }
    }

    if (got_pair_support) {
      ++m_e->m_key.pair_match_count;
    }

    CHECK(m_path.longest_read_id()) << m_e->describe(m_br);

    if (got_pair_support) {
      m_got_pair_support_or_mateless = true;
    } else if (!m_opts.readmap->has_mate(*m_path.longest_read_id())) {
      // Mateless.
      m_got_pair_support_or_mateless = true;
    }
  }

  void check_local_rejoin() {
    if ((m_path.size() - m_path.anchor_len()) <= m_opts.bidir_min_local_ref_overlap) {
      return;
    }

    aoffset_t outer_left_offset = m_br->right_push_view_offset();
    outer_left_offset -= (m_path.size() - m_path.anchor_len());
    auto ext = m_view->get_scaffold().split_extent_at(outer_left_offset - 1);
    if (ext.second.size() < m_opts.bidir_min_local_ref_overlap) {
      return;
    }

    for (int adjust : {
             0,   // same size
             -1,  // single base insert
             1    // single base delete
         }) {
      dna_slice ref_slice = ext.second.subseq(adjust + 1, ext.second.size() - (adjust + 1));
      unsigned shared = ref_slice.shared_prefix_length(m_path.seq());
      if (shared > m_opts.bidir_min_local_ref_overlap) {
        if (m_br->try_rejoin(outer_left_offset + adjust, dna_slice(), m_path,
                             m_e->pair_match_count())) {
          if (m_trace) {
            std::cout << "Local rejoin successful, adjust=" << adjust << ", shared = " << shared
                      << "\n";
          }
          m_e->m_did_rejoin = true;
          return;
        }
      }
    }
  }

  void check_rejoin() {
    if (!m_cur_read_r.valid()) {
      return;
    }

    if (m_trace) {
      std::cout << "Checking for rejoin at " << m_path << "\n";
    }

    auto rc_ri_entries = m_rc_view->range_info().entries_starting_with(m_cur_read_rc_r);

    for (const auto& r_and_ri : rc_ri_entries) {
      const auto& rc_r = r_and_ri.first;
      const auto& rc_ri = r_and_ri.second;
      CHECK_GE(rc_r.begin(), m_cur_read_rc_r.begin());
      CHECK_LE(rc_r.end(), m_cur_read_rc_r.end());

      if (m_trace) {
        std::cout << "Candidate range info: " << rc_r.sequence().rev_comp() << ": " << rc_ri
                  << "\n";
      }

      for (aoffset_t rc_ref_loc : rc_ri.reference_offsets) {
        aoffset_t left_offset = m_view->reverse_offset(rc_ref_loc + m_cur_read_rc_r.size());

        if (m_br->try_rejoin(left_offset, dna_slice(), m_path, m_e->pair_match_count())) {
          if (m_trace) {
            std::cout << "Ref rejoin success at: " << left_offset << "\n";
          }
          m_e->m_did_rejoin = true;
        } else {
          if (m_trace) {
            std::cout << "Ref rejoin failure at: " << left_offset << "\n";
          }
          if (rc_r.size() == m_cur_read_rc_r.size()) {
            if (m_trace) {
              std::cout << "Unable to rejoin at full entry; marking m_ambiguous since this might "
                           "start tracing a vaguely homologous section of reference.\n";
            }
            m_ambig = true;
          }
          continue;
        }
      }

      for (const auto& rc_rp : rc_ri.right_partials) {
        aoffset_t outer_left_offset = m_view->reverse_offset(rc_rp.outer_right_offset);
        dna_slice rp_seq = dna_slice(rc_rp.seq).rev_comp();
        CHECK_GE(rp_seq.size(), m_cur_read_rc_r.size());
        // Remove the part of it we've already traced
        CHECK_EQ(rp_seq.subseq(rp_seq.size() - m_cur_read_rc_r.size(), m_cur_read_rc_r.size()),
                 m_path.seq().subseq(0, m_cur_read_rc_r.size()))
            << m_e->describe(m_br);
        rp_seq = rp_seq.subseq(0, rp_seq.size() - m_cur_read_rc_r.size());
        if (m_br->try_rejoin(outer_left_offset, rp_seq, m_path, m_e->pair_match_count())) {
          if (m_trace) {
            std::cout << "Successfully rejoined right partial: " << rc_rp << "\n";
          }
          m_e->m_did_rejoin = true;
        } else {
          if (m_trace) {
            std::cout << "Failed to rejoin right partial: " << rc_rp << "\n";
          }
        }
      }
    }

    if (m_view->opts().rmap && m_view->opts().rmap->get(m_path.range().begin()).match_count() > 0) {
      // We're matching reference somewhere, not just somewhere that's valid to
      // rejoin.  Make sure we don't assemble the reference from
      // elsewhere when we're trying to trace this.
      m_ambig = true;
      if (m_trace) {
        std::cout << "Found ref without rejoin; marking as m_ambiguous.\n";
      }
    }
  }

  bool trace_one_base() WARN_UNUSED {
    size_t best_size = 0;
    dna_base best_base;
    seqset_range best_pushed;

    int best_size_count = 0;

    for (dna_base b : dna_bases()) {
      seqset_range pushed =
          m_path.range().push_front_drop(b, m_opts.min_overlap + m_path.bases_since_read());
      if (pushed.valid()) {
        if (pushed.size() == best_size) {
          ++best_size_count;
        } else if (pushed.size() > best_size) {
          best_size_count = 1;
          best_size = pushed.size();
          best_base = b;
          best_pushed = pushed;
        }
      }
    }

    if (best_size_count > 1) {
      if (m_trace) {
        std::cout << best_size_count << " paths found forward at pushed size " << best_pushed.size()
                  << "; setting ambiguous\n";
      }
      m_ambig = true;
      return false;
    }

    if (best_size_count == 0) {
      if (m_trace) {
        std::cout << best_size_count << " paths found forward; setting ambiguous\n";
      }
      m_ambig = true;
      return false;
    }

    CHECK_EQ(m_path.path_overlap(), m_e->m_key.path_overlap);
    m_path.push_front_drop(best_base, best_pushed);
    CHECK_LE(m_path.path_overlap(), m_e->m_key.path_overlap);
    m_e->m_key.path_overlap = m_path.path_overlap();
    if (m_path.loop_detected()) {
      if (m_trace) {
        std::cout << "Loop detected in path; setting ambiguous: " << m_path << "\n";
      }
      m_ambig = true;
      return false;
    }
    m_e->m_did_rejoin = false;
    return true;
  }

  push_search_entry* const m_e;
  branch* const m_br;
  state* const m_st;
  view_t* const m_view;
  view_t* const m_rc_view;
  const assemble_options& m_opts;
  path& m_path;

  seqset_range m_cur_read_r;
  seqset_range m_cur_read_rc_r;

  std::unique_ptr<pop_search_entry> m_pending_pop;
  boost::optional<right_partial> m_pending_rp;
  boost::optional<seqset_range> m_pending_rp_r;

  std::unordered_set<seqset_range, seqset_range_hash> m_seen;

  bool m_ambig = false;
  bool m_trace = k_trace_all;
  bool m_got_pair_support_or_mateless = false;
};

search_result push_search_entry::search_internal(branch* br) {
  push_tracer tr(this, br);
  return tr.search();
}

void push_search_entry::check_invariants(const branch* br) const {
  br->check_path_invariants(m_path);
}

std::string push_search_entry::describe_internal(const branch* br) const {
  std::stringstream result;

  view_t* v = br->push_view();
  if (v->is_rev_comp()) {
    result << "rev-push: " << br->opts().scaffold_name << ":"
           << v->reverse_offset(br->right_push_view_offset()) << " -> ";
    m_path.display_rev(result);
    result << "\n";
  } else {
    result << "fwd-push: ";
    m_path.display_fwd(result);
    result << " -> " << br->opts().scaffold_name << ":" << br->right_push_view_offset();
    result << "\n";
  }

  return result.str();
}

push_search_entry::push_search_entry(path p, unsigned pair_match_count)
    : branch_search_entry(search_entry_key(search_priority::PUSH, p, pair_match_count)) {
  m_path = std::move(p);
}

}  // namespace discovery
}  // namespace variants

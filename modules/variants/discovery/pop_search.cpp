#include "modules/variants/discovery/pop_search.h"

#include "modules/variants/discovery/push_search.h"
#include "modules/variants/discovery/rejoin.h"

namespace variants {
namespace discovery {
namespace {

constexpr bool k_trace_all = false;
constexpr bool k_trace_selected = true;

}  // namespace

class pop_tracer {
 public:
  pop_tracer(pop_search_entry* e, branch* br)
      : m_e(e),
        m_br(br),
        m_view(m_br->pop_view()),
        m_opts(m_view->opts()),
        m_left_offset(m_br->left_pop_view_offset()),
        m_max_right_offset(m_br->pop_view_farthest_right_offset()),
        m_popped(m_e->m_popped),
        m_rc_path(m_e->m_rc_path) {}

  search_result search() {
    if (k_trace_selected && m_br->trace_enabled(m_rc_path)) {
      m_trace = true;
    }

    if (m_rc_path.range().size() != m_popped.size()) {
      if (!m_br->explore(m_popped)) {
        return search_result::STOP_SEARCHING;
      }
    }

    if (k_trace_all) {
      std::cout << "POP search looking up: " << m_popped.sequence() << "\n";
      std::cout << "Path overlap: " << m_e->m_key.path_overlap << "\n";
    }

    bool found_any = false;
    for (const auto& r_and_ri : m_view->range_info().entries_starting_with(m_popped)) {
      const auto& r = r_and_ri.first;
      const auto& ri = r_and_ri.second;

      found_any = true;
      if (m_trace) {
        std::cout << "Considering range " << r.sequence() << ": " << ri << "\n";
      }

      for (aoffset_t ref_offset : ri.reference_offsets) {
        aoffset_t outer_ref_offset = ref_offset + m_popped.size();
        aoffset_t rc_outer_left_offset = m_view->reverse_offset(outer_ref_offset);

        if (m_trace) {
          aoffset_t seq_size = m_rc_path.size();
          aoffset_t ref_size = outer_ref_offset - m_left_offset;
          std::cout << "Considering rejoining ref at ref offset = " << ref_offset
                    << ", outer ref offset = " << outer_ref_offset
                    << ", rc outer left offset = " << rc_outer_left_offset
                    << " refsize = " << ref_size << " seqsize = " << seq_size
                    << " svlen=" << seq_size - ref_size << "\n";
        }

        m_br->try_rejoin(rc_outer_left_offset, dna_slice(), m_rc_path, m_e->pair_match_count());
      }

      if (r.size() > m_popped.size()) {
        interval_t relevant(m_left_offset, m_max_right_offset);
        interval_set_t relevant_supported_offsets = ri.pair_supported_offsets & relevant;
        if (!relevant_supported_offsets.empty()) {
          if (m_trace) {
            std::cout << "Making new search for pair-supported offsets "
                      << relevant_supported_offsets << "\n";
          }
          if (m_br->opts().bidir_pop_makes_push) {
            make_push_search_from_pair(r);
          } else {
            make_pop_search_from_pair(r);
          }
        } else if (!ri.pair_supported_offsets.empty()) {
          if (m_trace) {
            std::cout << "Pair support, but too far from our relevant offsets " << relevant << "\n";
          }
        }
      } else {
        if (m_trace) {
          std::cout << "Found range info at current range; not checking pair support";
        }
      }

      for (const auto& rp : ri.right_partials) {
        aoffset_t rc_outer_left_offset = m_view->reverse_offset(rp.outer_right_offset);
        dna_slice rp_seq = rp.seq;
        if (m_opts.bidir_validate_trace_state) {
          CHECK_EQ(rp_seq.subseq(0, m_popped.size()), m_popped.sequence());
        }
        // Trim off the bases that are repeated in m_popped
        rp_seq = rp_seq.subseq(m_popped.size(), rp.seq.size() - m_popped.size());
        if (m_trace) {
          std::cout << "Considering rejoining right partial: " << rp << "\n";
        }
        m_br->try_rejoin(rc_outer_left_offset, rp_seq.rev_comp(), m_rc_path,
                         m_e->pair_match_count() + rp.pair_match_count);
      }
    }

    if (found_any && m_trace) {
      std::cout << "Done considering all pop rejoin ranges.\n";
    }

    uint64_t popped_diff = m_popped.end() - m_popped.begin();
    uint64_t seqset_size = m_opts.seqset->size();
    bool popped_too_much_of_seqset = popped_diff * m_opts.bidir_max_pop_seqset_portion > seqset_size;
    if (popped_too_much_of_seqset && m_trace) {
      std::cout << "Popped is too general general: " << popped_diff << " vs seqset size "
                << seqset_size << "; ratio is 1:" << (double(seqset_size) / popped_diff) << " ("
                << printstring("%.6f%%", popped_diff * 100. / seqset_size) << ")\n"
                << "Seq: " << m_popped.sequence() << "\n";
    }
    if (m_popped.size() > m_opts.min_pop_overlap && !popped_too_much_of_seqset) {
      m_popped = m_popped.pop_front();
      m_e->m_key.path_overlap = std::min<unsigned>(m_e->m_key.path_overlap, m_popped.size());
      m_e->m_key.tot_overlap += m_popped.size();
      ++m_e->m_key.tot_overlap_bases;
      if (k_trace_all) {
        std::cout << "New popped: " << m_popped.sequence() << "\n";
        std::cout << "New path overlap: " << m_e->m_key.path_overlap << "\n";
      }
      return search_result::SEARCH_MORE;
    } else {
      if (m_trace) {
        std::cout << "Popped is " << m_popped.sequence() << " of size " << m_popped.size()
                  << ", too small to pop more\n";
      }
      return search_result::STOP_SEARCHING;
    }
  }

  void make_pop_search_from_pair(const seqset_range& r) {
    CHECK_GE(r.begin(), m_popped.begin());
    CHECK_LE(r.end(), m_popped.end());
    CHECK_GT(r.size(), m_popped.size());

    auto maybe_rd = m_opts.readmap->get_longest_prefix_read(r);
    CHECK(maybe_rd) << r.sequence();
    readmap::read rd = *maybe_rd;

    CHECK_EQ(rd.size(), r.size());

    seqset_range rc_r = rd.get_rev_comp().get_seqset_entry();
    CHECK_EQ(r.size(), rc_r.size());

    path rc_push_path = m_rc_path;

    rc_push_path.push_front_drop(rc_r.sequence(r.size() - m_popped.size()));

    if (rc_push_path.loop_detected()) {
      return;
    }

    // Save right partial
    aoffset_t rc_outer_right_offset = m_br->right_push_view_offset() + rc_push_path.anchor_len();
    right_partial rc_rp(rc_push_path.seq(), rc_outer_right_offset, m_e->pair_match_count() + 1);
    m_view->reverse_view()->add_right_partial(rc_r, std::move(rc_rp));

    std::unique_ptr<pop_search_entry> pop =
        make_unique<pop_search_entry>(r, std::move(rc_push_path), m_e->pair_match_count() + 1);
    if (m_opts.bidir_validate_trace_state) {
      pop->check_invariants(m_br);
    }
    m_br->add_search_entry(std::move(pop));
  }

  void make_push_search_from_pair(const seqset_range& r) {
    CHECK_GE(r.begin(), m_popped.begin());
    CHECK_LE(r.end(), m_popped.end());
    CHECK_GT(r.size(), m_popped.size());

    auto maybe_rd = m_opts.readmap->get_longest_prefix_read(r);
    CHECK(maybe_rd) << r.sequence();
    readmap::read rd = *maybe_rd;

    CHECK_EQ(rd.size(), r.size());

    seqset_range rc_r = rd.get_rev_comp().get_seqset_entry();
    CHECK_EQ(r.size(), rc_r.size());

    path rc_push_path = m_rc_path;

    rc_push_path.push_front_drop(rc_r.sequence(r.size() - m_popped.size()));

    if (rc_push_path.loop_detected()) {
      return;
    }
    std::unique_ptr<push_search_entry> push =
        make_unique<push_search_entry>(std::move(rc_push_path), m_e->pair_match_count() + 1);
    if (m_opts.bidir_validate_trace_state) {
      push->check_invariants(m_br);
    }
    m_br->add_search_entry(std::move(push));
  }

 private:
  pop_search_entry* const m_e;
  branch* const m_br;
  view_t* const m_view;
  const assemble_options& m_opts;
  aoffset_t const m_left_offset;
  aoffset_t const m_max_right_offset;
  seqset_range& m_popped;
  const path& m_rc_path;
  bool m_trace = k_trace_all;
};

search_result pop_search_entry::search_internal(branch* br) {
  pop_tracer tr(this, br);
  return tr.search();
}

void pop_search_entry::check_invariants(const branch* br) const {
  CHECK_GE(m_rc_path.size(), m_popped.size()) << describe(br);
  dna_slice seq = m_rc_path.seq().rev_comp();
  if (br->pop_view()->opts().bidir_validate_trace_state > 1) {
    CHECK_EQ(seq.subseq(seq.size() - m_popped.size(), m_popped.size()), m_popped.sequence())
        << "\nSeq: " << seq << "\nR: " << m_popped.sequence() << "\nSearch entry: " << describe(br);
  }

  br->check_path_invariants(m_rc_path);
}

std::string pop_search_entry::describe_internal(const branch* br) const {
  std::stringstream result;

  view_t* v = br->pop_view();
  if (v->is_rev_comp()) {
    result << "rev-POP@" << v->reverse_offset(br->left_pop_view_offset());
  } else {
    result << "fwd-POP@" << br->left_pop_view_offset();
  }

  result << ":  " << m_popped.sequence() << " for reverse path " << m_rc_path << "\n";
  return result.str();
}

pop_search_entry::pop_search_entry(const seqset_range& popped, path rc_path,
                                   unsigned pair_match_count)
    : branch_search_entry(search_entry_key(search_priority::POP, rc_path, pair_match_count)),
      m_popped(popped) {
  m_rc_path = std::move(rc_path);
}

}  // namespace discovery
}  // namespace variants

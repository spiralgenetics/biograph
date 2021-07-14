#pragma once

#include "modules/variants/assemble_testutil.h"
#include "modules/variants/discovery/pop_search.h"
#include "modules/variants/discovery/push_search.h"
#include "modules/variants/discovery/rejoin.h"
#include "modules/variants/discovery/state.h"
#include "modules/variants/discovery/view.h"
#include "modules/variants/discovery/walk_ref.h"

namespace variants {
namespace discovery {

class discovery_test : public assemble_test {
 public:
  // Data access for search entries, since we're a friend class.
  template<typename T>
  unsigned path_overlap(const T* e) const {
    return e->get_key().m_path_overlap;
  }

  bool is_push_rev_comp(const branch* br) const { return br->push_view()->is_rev_comp(); }
  bool is_pop_rev_comp(const branch* br) const { return br->pop_view()->is_rev_comp(); }
  seqset_range range(const push_search_entry* e) const { return e->m_path.range(); }
  dna_slice seq(const push_search_entry* e) const { return e->m_path.seq(); }
  seqset_range range(const pop_search_entry* e) const { return e->m_popped; }
  dna_slice seq(const pop_search_entry* e) const {
    return dna_slice(e->m_rc_path.seq()).rev_comp();
  }

  aoffset_t left_offset(const rejoin_search_entry* e) const { return e->m_left_offset; }
  dna_slice seq(const rejoin_search_entry* e) const { return e->m_path.seq(); }

  const range_info_table_t& range_info_table(view_t* v) const { return v->m_range_info; }

  std::vector<const branch_search_entry*> search_entries(const branch* br) const {
    std::vector<const branch_search_entry*> result;
    for (const auto& e : br->m_search_entries) {
      result.push_back(e.get());
    }
    return result;
  }

  void execute_search(branch* br, branch_search_entry_ptr e) {
    br->execute_search_for_testing(std::move(e));
    m_st->check_invariants();
  }
  bool search_each_branch_once() {
    bool did_any = false;
    for (view_t* v : {fwd_view(), rev_view()}) {
      for (branch* br : branches(v)) {
        if (!br->empty()) {
          br->execute_one_search_for_testing();
          did_any = true;
        }
      }
    }
    m_st->check_invariants();
    return did_any;
  }

  std::vector<branch*> branches(view_t* v) const {
    std::vector<branch*> result;
    for (const auto& br_item : v->m_branches) {
      result.push_back(br_item.second.get());
    }
    return result;
  }

// Have to reify all the branch-and-search-entry pairs so that gtest PrintTo overloads will work
// right. :(
#define REIFY_PAIR(PAIR_TYPE, ENTRY_TYPE)                                                          \
  struct PAIR_TYPE {                                                                               \
    const branch* first;                                                                           \
    const ENTRY_TYPE* second;                                                                      \
    PAIR_TYPE() = default;                                                                         \
    PAIR_TYPE(const branch* new_first, const ENTRY_TYPE* new_second)                               \
        : first(new_first), second(new_second) {}                                                  \
    PAIR_TYPE(std::pair<const branch*, const ENTRY_TYPE*> p) : first(p.first), second(p.second) {} \
  };
  REIFY_PAIR(push_entry_pair, push_search_entry);
  REIFY_PAIR(pop_entry_pair, pop_search_entry);
  REIFY_PAIR(rejoin_entry_pair, rejoin_search_entry);
#undef REIFY_PAIR

  view_t* fwd_view() { return m_rev_comp ? m_st->rev_view() : m_st->fwd_view(); }
  view_t* rev_view() { return fwd_view()->reverse_view(); }

 protected:
  ~discovery_test() {
    if (m_st) {
      m_st->check_invariants();
    }
  }
  void init_discovery() {
    m_fwd_scaffold = *m_options.scaffold;
    m_rev_scaffold = m_fwd_scaffold.rev_comp();
    if (m_rev_comp) {
      m_options.scaffold = &m_rev_scaffold;
    }
    m_options.bidir_validate_trace_state = 1000;
    m_st.emplace(m_options, test_output());
  }

  ::testing::Matcher<push_entry_pair> FwdPushSearchEntry(unsigned path_overlap,
                                                         aoffset_t right_offset, dna_sequence seq,
                                                         dna_sequence r);
  ::testing::Matcher<push_entry_pair> RevPushSearchEntry(unsigned path_overlap,
                                                         aoffset_t left_offset, dna_sequence seq,
                                                         dna_sequence r);
  ::testing::Matcher<push_entry_pair> PushSearchEntry(bool is_rev_comp, unsigned path_overlap,
                                                      aoffset_t right_offset, dna_sequence seq,
                                                      dna_sequence range);
  ::testing::Matcher<rejoin_entry_pair> RejoinSearchEntry(unsigned path_overlap,
                                                          aoffset_t left_offset, dna_sequence seq,
                                                          aoffset_t right_offset);
  ::testing::Matcher<pop_entry_pair> FwdPopSearchEntry(unsigned path_overlap, aoffset_t left_offset,
                                                       const dna_sequence& seq,
                                                       const dna_sequence& r);
  ::testing::Matcher<pop_entry_pair> RevPopSearchEntry(unsigned path_overlap,
                                                       aoffset_t right_offset,
                                                       const dna_sequence& seq,
                                                       const dna_sequence& r);
  ::testing::Matcher<pop_entry_pair> PopSearchEntry(bool is_rev_comp, unsigned path_overlap,
                                                    aoffset_t left_offset, const dna_sequence& seq,
                                                    const dna_sequence& r);

  seqset_range get_seqset_range(dna_slice seq) {
    seqset_range r = m_st->m_options.seqset->find(seq);
    CHECK(r.valid()) << seq;
    return r;
  }

  void save_partials() {
    CHECK(m_st);
    m_st->check_invariants();
    m_left_partials.clear();
    m_right_partials.clear();

    for (view_t* v : m_st->both_dirs()) {
      bool is_rev_comp = v->is_rev_comp() != m_rev_comp;

      auto& partials = is_rev_comp ? m_left_partials : m_right_partials;
      for (const auto& r_and_ri : v->range_info()) {
        const auto& ri = r_and_ri.second;

        for (const auto& rp : ri.right_partials) {
          partials.emplace(
              is_rev_comp ? v->reverse_offset(rp.outer_right_offset) : (rp.outer_right_offset),
              is_rev_comp ? rp.seq.rev_comp() : rp.seq);
        }
      }
    }
  }

  void save_search_entries();

  void save_pair_support() {
    CHECK(m_st);
    m_st->check_invariants();

    m_pair_support.clear();
    m_rev_pair_support.clear();

    for (view_t* v : m_st->both_dirs()) {
      bool is_rev_comp = v->is_rev_comp() != m_rev_comp;
      auto& pair_support_table = is_rev_comp ? m_rev_pair_support : m_pair_support;
      for (const auto& r_and_ri : range_info_table(v)) {
        const seqset_range& r = r_and_ri.first;
        const range_info_t& ri = r_and_ri.second;

        if (ri.pair_supported_offsets.empty()) {
          continue;
        }
        dna_sequence seq = r.sequence();

        interval_set_t offsets = ri.pair_supported_offsets;
        if (is_rev_comp) {
          seq = seq.rev_comp();
          interval_set_t rev_offsets;
          for (const auto& offset : offsets) {
            rev_offsets += interval_t(v->reverse_offset(offset.upper()) - seq.size(),
                                      v->reverse_offset(offset.lower()) - seq.size());
          }
          offsets = rev_offsets;
        }
        pair_support_table[seq] += offsets;
      }
    }

    EXPECT_THAT(m_pair_support, ::testing::ContainerEq(m_rev_pair_support));
  }

  void add_ref() {
    for (view_t* v : m_st->both_dirs()) {
      walk_ref_t wr(v);
      wr.walk_ref(0, v->get_scaffold().end_pos());
      wr.check_invariants();
      m_st->check_invariants();
      wr.init_pairs_and_push();
      m_st->check_invariants();
    }
  }

  void add_ref_without_search() {
    for (view_t* v : m_st->both_dirs()) {
      walk_ref_t wr(v);
      wr.walk_ref(0, v->get_scaffold().end_pos());
      wr.check_invariants();
      m_st->check_invariants();
      wr.init_pairs_and_push();
      m_st->check_invariants();
      m_st->discard_search_entries();
    }
  }

 protected:
  boost::optional<state> m_st;
  scaffold m_fwd_scaffold;
  scaffold m_rev_scaffold;
  bool m_rev_comp = false;

  std::vector<push_entry_pair> m_push_entries;
  std::vector<pop_entry_pair> m_pop_entries;
  std::vector<rejoin_entry_pair> m_rejoin_entries;

  std::map<dna_sequence, interval_set_t> m_pair_support;
  std::map<dna_sequence, interval_set_t> m_rev_pair_support;

  std::set<std::pair<aoffset_t, dna_sequence>> m_left_partials;
  std::set<std::pair<aoffset_t, dna_sequence>> m_right_partials;
};

#define DEFINE_PRINT(PAIR_TYPE)                               \
  inline void PrintTo(PAIR_TYPE br_and_e, std::ostream* os) { \
    (*os) << br_and_e.second->describe(br_and_e.first);       \
  }
DEFINE_PRINT(const discovery_test::push_entry_pair&);
DEFINE_PRINT(const discovery_test::pop_entry_pair&);
DEFINE_PRINT(const discovery_test::rejoin_entry_pair&);
#undef DEFINE_PRINT

inline void PrintTo(const view_t* v, std::ostream* os) {
  if (v->is_rev_comp()) {
    (*os) << "View(rev)";
  } else {
    (*os) << "View(fwd)";
  }
}
inline void PrintTo(view_t* v, std::ostream* os) { (*os) << (const view_t*)v; }

}  // namespace discovery
}  // namespace variants

// TODO(nils): Is there a way to resolve to these printers besides putting this in the boost::icl
// namespace?
namespace boost {
namespace icl {

inline void PrintTo(const variants::discovery::interval_t& intv, std::ostream* os) {
  (*os) << intv;
}
inline void PrintTo(const variants::discovery::interval_set_t& intvs, std::ostream* os) {
  (*os) << intvs;
}
inline void PrintTo(const variants::discovery::ploids_remaining_t& remaining, std::ostream* os) {
  (*os) << remaining;
}

}  // namespace icl
}  // namespace boost

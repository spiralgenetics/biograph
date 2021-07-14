#include "modules/variants/discovery/discovery_testutil.h"

using namespace testing;

namespace variants {
namespace discovery {

namespace {

// Canonicalize orders for search entries so we can compare with
// ElementsAre instead of UnorderedElementsAre.
template <typename T>
struct canon_order {};

template <typename PairT, typename T>
struct br_canon_order {
  br_canon_order(discovery_test* t) : m_t(t) {}

  bool operator()(const PairT& a, const PairT& b) const {
    if (a.first->push_view()->is_rev_comp() != b.first->push_view()->is_rev_comp()) {
      return a.first->push_view()->is_rev_comp() == m_t->fwd_view()->is_rev_comp();
    }

    if (a.first->right_push_view_offset() != b.first->right_push_view_offset()) {
      return a.first->right_push_view_offset() < b.first->right_push_view_offset();
    }

    if (a.second->path_overlap() != b.second->path_overlap()) {
      return a.second->path_overlap() > b.second->path_overlap();
    }

    return canon_order<T>()(m_t, a.second, b.second);
  }

  discovery_test* m_t = nullptr;
};

template <>
struct canon_order<push_search_entry> {
  bool operator()(discovery_test* t, const push_search_entry* a, const push_search_entry* b) {
    return t->seq(a) < t->seq(b);
  }
};

template <>
struct canon_order<pop_search_entry> {
  bool operator()(discovery_test* t, const pop_search_entry* a, const pop_search_entry* b) {
    if (t->seq(a) != t->seq(b)) {
      return t->seq(a) < t->seq(b);
    }
    return t->range(a).size() < t->range(b).size();
  }
};

template <>
struct canon_order<rejoin_search_entry> {
  bool operator()(discovery_test* t, const rejoin_search_entry* a, const rejoin_search_entry* b) {
    if (t->left_offset(a) != t->left_offset(b)) {
      return t->left_offset(a) < t->left_offset(b);
    }
    return t->seq(a) < t->seq(b);
  }
};

}  // namespace

void discovery_test::save_search_entries() {
  CHECK(m_st);
  m_st->check_invariants();

  m_push_entries.clear();
  m_pop_entries.clear();
  m_rejoin_entries.clear();

  for (bool rc : {false, true}) {
    view_t* v = rc ? fwd_view() : rev_view();
    for (const auto& br_item : v->m_branches) {
      const branch* br = br_item.second.get();
      for (const branch_search_entry* e : search_entries(br)) {
        e->check_invariants(br);

        const auto* push = dynamic_cast<const push_search_entry*>(e);
        if (push) {
          m_push_entries.push_back(std::make_pair(br, push));
        }

        const auto* pop = dynamic_cast<const pop_search_entry*>(e);
        if (pop) {
          m_pop_entries.push_back(std::make_pair(br, pop));
        }

        const auto* rejoin = dynamic_cast<const rejoin_search_entry*>(e);
        if (rejoin) {
          m_rejoin_entries.push_back(std::make_pair(br, rejoin));
        }

        CHECK(push || pop || rejoin) << "Unknown search entry type: " << e->describe(br);
      }
    }
  }

  // Put them in a canonical order so tests don't have to use Unordered.
  std::sort(m_push_entries.begin(), m_push_entries.end(),
            br_canon_order<push_entry_pair, push_search_entry>(this));
  std::sort(m_pop_entries.begin(), m_pop_entries.end(),
            br_canon_order<pop_entry_pair, pop_search_entry>(this));
  std::sort(m_rejoin_entries.begin(), m_rejoin_entries.end(),
            br_canon_order<rejoin_entry_pair, rejoin_search_entry>(this));
}

MATCHER_P6(PushSearchEntryInternal, t, is_rev_comp, path_overlap, right_offset, seq, range, "") {
  const branch* br = arg.first;
  const push_search_entry* e = dynamic_cast<const push_search_entry*>(arg.second);

  if (!e) {
    *result_listener << " where the search entry is not a push search entry";
    return false;
  }

  bool ok = true;

  if (e->path_overlap() != path_overlap) {
    ok = false;
    *result_listener << " where path_overlap " << e->path_overlap() << " is not expected "
                     << path_overlap << "\n";
  }
  if (t->is_push_rev_comp(br) != is_rev_comp) {
    ok = false;
    *result_listener << " where is_rev_comp " << t->is_push_rev_comp(br) << " is not expected "
                     << is_rev_comp << "\n";
  }
  if (aoffset_t(right_offset) != br->right_push_view_offset()) {
    ok = false;
    *result_listener << " where offset " << br->right_push_view_offset() << " is not expected "
                     << right_offset << "\n";
  }
  if (range != t->range(e).sequence()) {
    ok = false;
    *result_listener << " where range " << t->range(e).sequence() << " is not expected " << range
                     << "\n";
  }
  if (seq != t->seq(e)) {
    ok = false;
    *result_listener << " where sequence " << t->seq(e) << " is not expected " << seq << "\n";
  }

  return ok;
}

::testing::Matcher<discovery_test::push_entry_pair> discovery_test::PushSearchEntry(
    bool is_rev_comp, unsigned path_overlap, aoffset_t right_offset, dna_sequence seq,
    dna_sequence range) {
  if (m_rev_comp) {
    is_rev_comp = !is_rev_comp;
  }

  return PushSearchEntryInternal(this, is_rev_comp, path_overlap, right_offset, seq, range);
}

::testing::Matcher<discovery_test::push_entry_pair> discovery_test::FwdPushSearchEntry(
    unsigned path_overlap, aoffset_t right_offset, dna_sequence seq, dna_sequence range) {
  return PushSearchEntry(false /* !rev_comp */, path_overlap, right_offset, seq, range);
}

::testing::Matcher<discovery_test::push_entry_pair> discovery_test::RevPushSearchEntry(
    unsigned path_overlap, aoffset_t left_offset, dna_sequence seq, dna_sequence range) {
  return PushSearchEntry(true /* rev_comp */, path_overlap,
                         m_st->fwd_view()->reverse_offset(left_offset), seq.rev_comp(),
                         range.rev_comp());
}

MATCHER_P6(RejoinSearchEntryInternal, t, is_rev_comp, path_overlap, left_offset, seq, right_offset,
           "") {
  const branch* br = arg.first;
  const rejoin_search_entry* e = dynamic_cast<const rejoin_search_entry*>(arg.second);

  if (!e) {
    *result_listener << " where the search entry is not a rejoin search entry";
    return false;
  }

  bool ok = true;

  if (e->path_overlap() != path_overlap) {
    ok = false;
    *result_listener << " where path_overlap " << e->path_overlap() << " is not expected "
                     << path_overlap << "\n";
  }
  if (t->is_push_rev_comp(br) != is_rev_comp) {
    ok = false;
    *result_listener << " where is_rev_comp " << t->is_push_rev_comp(br) << " is not expected "
                     << is_rev_comp << "\n";
  }
  if (aoffset_t(left_offset) != t->left_offset(e)) {
    ok = false;
    *result_listener << " where left offset " << t->left_offset(e) << " is not expected "
                     << left_offset << "\n";
  }
  if (seq != t->seq(e)) {
    ok = false;
    *result_listener << " where sequence " << t->seq(e) << " is not expected " << seq << "\n";
  }
  if (aoffset_t(right_offset) != br->right_push_view_offset()) {
    ok = false;
    *result_listener << " where right offset " << br->right_push_view_offset()
                     << " is not expected " << right_offset << "\n";
  }

  return ok;
}

::testing::Matcher<discovery_test::rejoin_entry_pair> discovery_test::RejoinSearchEntry(
    unsigned path_overlap, aoffset_t left_offset, dna_sequence seq, aoffset_t right_offset) {
  return AnyOf(
      RejoinSearchEntryInternal(this, m_rev_comp, path_overlap, left_offset, seq, right_offset),
      RejoinSearchEntryInternal(this, !m_rev_comp, path_overlap,
                                m_st->rev_view()->reverse_offset(right_offset), seq.rev_comp(),
                                m_st->rev_view()->reverse_offset(left_offset)));
}

MATCHER_P6(PopSearchEntryInternal, t, is_rev_comp, path_overlap, left_offset, seq, range, "") {
  const branch* br = arg.first;
  const pop_search_entry* e = dynamic_cast<const pop_search_entry*>(arg.second);

  if (!e) {
    *result_listener << " where the search entry is not a pop search entry";
    return false;
  }

  bool ok = true;

  if (e->path_overlap() != path_overlap) {
    ok = false;
    *result_listener << " where path_overlap " << e->path_overlap() << " is not expected "
                     << path_overlap << "\n";
  }
  if (t->is_pop_rev_comp(br) != is_rev_comp) {
    ok = false;
    *result_listener << " where is_rev_comp " << t->is_pop_rev_comp(br) << " is not expected "
                     << is_rev_comp << "\n";
  }
  if (aoffset_t(left_offset) != br->left_pop_view_offset()) {
    ok = false;
    *result_listener << " where left offset " << br->left_pop_view_offset() << " is not expected "
                     << left_offset << "\n";
  }

  if (range != t->range(e).sequence()) {
    ok = false;
    *result_listener << " where range " << t->range(e).sequence() << " is not expected " << range
                     << "\n";
  }
  if (seq != t->seq(e)) {
    ok = false;
    *result_listener << " where sequence " << t->seq(e) << " is not expected " << seq << "\n";
  }

  return ok;
}

::testing::Matcher<discovery_test::pop_entry_pair> discovery_test::PopSearchEntry(
    bool is_rev_comp, unsigned path_overlap, aoffset_t left_offset, const dna_sequence& seq,
    const dna_sequence& range) {
  if (m_rev_comp) {
    is_rev_comp = !is_rev_comp;
  }
  return PopSearchEntryInternal(this, is_rev_comp, path_overlap, left_offset, seq, range);
}

::testing::Matcher<discovery_test::pop_entry_pair> discovery_test::FwdPopSearchEntry(
    unsigned path_overlap, aoffset_t left_offset, const dna_sequence& seq,
    const dna_sequence& range) {
  return PopSearchEntry(false /* !rev_comp */, path_overlap, left_offset, seq, range);
}

::testing::Matcher<discovery_test::pop_entry_pair> discovery_test::RevPopSearchEntry(
    unsigned path_overlap, aoffset_t right_offset, const dna_sequence& seq,
    const dna_sequence& range) {
  return PopSearchEntry(true /* rev_comp */, path_overlap,
                        m_st->fwd_view()->reverse_offset(right_offset), seq.rev_comp(),
                        range.rev_comp());
}

}  // namespace discovery
}  // namespace variants

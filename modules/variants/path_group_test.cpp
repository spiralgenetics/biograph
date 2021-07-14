#include "modules/variants/path_group.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class path_group_test : public assemble_test, public path_group::listener {
 public:
  class mock_dobj : public path_group::distant_object {
   public:
    mock_dobj(path_group_test* t) : m_test(t) { ++m_test->m_dobj_created; }
    ~mock_dobj() override {
      EXPECT_LT(m_test->m_dobj_destroyed, m_test->m_dobj_created);
      ++m_test->m_dobj_destroyed;
    }

   private:
    path_group_test* m_test = nullptr;
  };

  class visitor : public path_group::dobj_visitor {
   public:
    visitor(path_group_test* t, const seqset_range& r, int* expected_visits)
        : m_test(t), m_r(r), m_expected_visits(expected_visits) {}
    void visit(path_group::distant_object*, int distance) override {
      std::cout << "Visiting distant object at " << distance << "\n";
      m_test->m_events.push_back({m_r, distance});

      EXPECT_GT(*m_expected_visits, 0);
      --*m_expected_visits;
    }

   private:
    path_group_test* m_test = nullptr;
    seqset_range m_r;
    int* m_expected_visits = nullptr;
  };

  void on_seqset_entry(const seqset_range& r, path_group* pg) override {
    if (!r.is_maximal()) {
      return;
    }
    std::cout << "On seqset entry: " << r.sequence() << "\n";
    visitor v(this, r, &m_expected_visits);
    pg->visit_distant_objects(r, v);
  }

  void TearDown() override {
    EXPECT_EQ(m_dobj_created, m_dobj_destroyed);
    assemble_test::TearDown();
  }

  struct event {
    seqset_range r;
    int distance;

    friend std::ostream& operator<<(std::ostream& os, const event& ev) {
      return os << "\nevent: r=" << ev.r.sequence().rev_comp() << " distance=" << ev.distance;
    }
  };
  std::vector<event> m_events;
  int m_dobj_created = 0;
  int m_dobj_destroyed = 0;
  constexpr static unsigned k_min_overlap = 20;
  int m_expected_visits = 1000;
};

constexpr unsigned path_group_test::k_min_overlap;

MATCHER_P2(EventIs, seq, distance, "") {
  return arg.r.sequence().rev_comp() == seq && int(arg.distance) == int(distance);
}

TEST_F(path_group_test, basic_listener) {
  use_reads({tseq("abcde"), tseq("cdefg")});
  auto pg = make_unique<path_group>(m_seqset->ctx_begin(), k_min_overlap, this);
  pg->add_distant_object(std::make_shared<mock_dobj>(this), 1);
  pg->add_sequence(tseq("abcdefg"));

  EXPECT_THAT(m_events, ElementsAre(EventIs(tseq("abcde"), tseq("abcde").size()),
                                    EventIs(tseq("cdefg"), tseq("abcdefg").size())));
}

TEST_F(path_group_test, basic_appends_in_parts) {
  use_reads({tseq("abcde"), tseq("cdefg")});
  auto pg = make_unique<path_group>(m_seqset->ctx_begin(), k_min_overlap, this);
  pg->add_sequence(tseq("ab"));
  pg->add_distant_object(std::make_shared<mock_dobj>(this), 1);
  pg->add_sequence(tseq("cdef"));
  pg->add_sequence(tseq("g"));

  EXPECT_THAT(m_events, ElementsAre(EventIs(tseq("abcde"), tseq("cde").size()),
                                    EventIs(tseq("cdefg"), tseq("cdefg").size())));
}

TEST_F(path_group_test, listener_expires) {
  use_reads({tseq("abcde") + dna_T, tseq("cde") + dna_T + tseq("fg")});
  auto pg =
      make_unique<path_group>(m_seqset->ctx_begin(), k_min_overlap, this);
  pg->add_distant_object(std::make_shared<mock_dobj>(this), tseq("abcde").size());
  pg->add_sequence(tseq("abcde"));
  EXPECT_THAT(m_events, IsEmpty());
  EXPECT_EQ(0, m_dobj_destroyed);
  pg->add_sequence(dna_T);
  EXPECT_THAT(m_events, ElementsAre(EventIs(tseq("abcde") + dna_T, tseq("abcde").size() + 1)));
  pg->flush();
  EXPECT_EQ(1, m_dobj_destroyed);
  pg->add_sequence(tseq("fg"));
  EXPECT_EQ(1, m_dobj_destroyed);
  pg.reset();
  EXPECT_EQ(1, m_dobj_destroyed);

  EXPECT_THAT(m_events, ElementsAre(EventIs(tseq("abcde") + dna_T, tseq("abcde").size() + 1)));
  EXPECT_EQ(1, m_dobj_destroyed);
}

TEST_F(path_group_test, paths_join) {
  use_reads({tseq("abcde"), tseq("cde") + dna_T + tseq("fg"), tseq("cdefg"), tseq("fghij")});
  auto pg = make_unique<path_group>(m_seqset->ctx_begin(), k_min_overlap, this);
  pg->add_distant_object(std::make_shared<mock_dobj>(this), 1);
  pg->add_sequence(tseq("abcde"));
  {
    auto split_off = pg->split();
    split_off->add_sequence(dna_T);
    pg->join(std::move(split_off));
  }
  EXPECT_EQ(2, pg->size());
  pg->add_sequence(tseq("fgh"));
  EXPECT_EQ(1, pg->size());
  pg->add_sequence(tseq("ij"));
  EXPECT_THAT(m_events, UnorderedElementsAre(EventIs(tseq("abcde"), tseq("abcde").size()),
                                             EventIs(tseq("cde") + dna_T + tseq("fg"),
                                                     tseq("abcde").size() + 1 + tseq("fg").size()),
                                             EventIs(tseq("cdefg"), tseq("abcdefg").size()),
                                             EventIs(tseq("fghij"), tseq("abcdefghij").size())));
}

TEST_F(path_group_test, DISABLED_min_overlap) {
  use_reads({tseq("abcde"), tseq("de") + dna_T + tseq("fgh"), tseq("e") + dna_A + tseq("fgh")});
  auto pg = make_unique<path_group>(m_seqset->ctx_begin(), k_min_overlap, this);
  pg->add_distant_object(std::make_shared<mock_dobj>(this), 1);
  {
    auto split_off = pg->split();
    pg->add_sequence(dna_T);
    split_off->add_sequence(dna_A);
    pg->join(std::move(split_off));
    EXPECT_EQ(2, pg->size());
  }
  pg->add_sequence(tseq("f"));
  EXPECT_EQ(1, pg->size());
  pg->add_sequence(tseq("gh"));
  pg.reset();
  EXPECT_THAT(m_events, ElementsAre(EventIs(tseq("abcde"), tseq("abcde").size()),
                                    EventIs(tseq("de") + dna_T + tseq("fgh"),
                                            tseq("abcde").size() + 1 + tseq("fgh").size())));
}

}  // namespace variants

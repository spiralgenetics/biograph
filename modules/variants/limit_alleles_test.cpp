#include "modules/variants/limit_alleles.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <functional>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class limit_alleles_test : public assemble_test {
 public:
  void add(aoffset_t left_offset, aoffset_t right_offset, dna_sequence seq, size_t prio) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = allocate_assembly_id();
    a->tags.insert("limit_alleles_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = seq;
    CHECK(!m_prios.count(seq));
    m_prios[seq] = prio;
    m_pipeline->add(std::move(a));

    CHECK(!m_prio_used.count(prio)) << "Priority should give a strict ordering: " << seq;
    m_prio_used.insert(prio);
  }

  void start(size_t max_alleles) {
    m_pipeline = make_unique<limit_alleles>(
        max_alleles, std::bind(&limit_alleles_test::sort_func, this, std::placeholders::_1),
        std::bind(&limit_alleles_test::limited_func, this, std::placeholders::_1), test_output());
  }

  void flush() {
    m_pipeline.reset();

    EXPECT_EQ(m_assemblies.size(), m_prios.size());
  }

 protected:
  std::vector<assembly_ptr> sort_func(std::vector<assembly_ptr> asms) {
    std::set<dna_sequence> sort_seqs;
    for (auto& a : asms) {
      CHECK(sort_seqs.insert(a->seq).second);
    }
    m_sorts.emplace_back(std::move(sort_seqs));

    sort(asms.begin(), asms.end(), [this](const assembly_ptr& a, const assembly_ptr& b) {
      CHECK(m_prios.count(a->seq)) << a->seq;
      CHECK(m_prios.count(b->seq)) << b->seq;
      return m_prios[a->seq] < m_prios[b->seq];
    });
    return asms;
  }

  void limited_func(const assembly_ptr& a) { m_limited.insert(a->seq); }

  std::unique_ptr<limit_alleles> m_pipeline;
  std::vector<std::set<dna_sequence>> m_sorts;
  std::set<dna_sequence> m_limited;
  std::map<dna_sequence, size_t> m_prios;
  std::set<size_t> m_prio_used;
};

TEST_F(limit_alleles_test, single_at_once) {
  start(1);
  add(1, 2, tseq("Single"), 1);
  flush();

  EXPECT_THAT(m_sorts, IsEmpty());
  EXPECT_THAT(m_limited, IsEmpty());
}

TEST_F(limit_alleles_test, simple_conflict) {
  start(1);
  add(1, 2, tseq("First"), 1);
  add(1, 2, tseq("Second"), 2);
  flush();

  EXPECT_THAT(m_sorts, ElementsAre(UnorderedElementsAre(tseq("First"), tseq("Second"))));
  EXPECT_THAT(m_limited, UnorderedElementsAre(tseq("Second")));
}

TEST_F(limit_alleles_test, two_at_once) {
  start(2);
  add(1, 2, tseq("First"), 1);
  add(1, 2, tseq("Second"), 2);
  flush();

  EXPECT_THAT(m_sorts, IsEmpty());
  EXPECT_THAT(m_limited, IsEmpty());
}

TEST_F(limit_alleles_test, two_not_at_once) {
  start(1);
  add(1, 2, tseq("First"), 1);
  add(2, 3, tseq("Second"), 2);
  flush();

  EXPECT_THAT(m_sorts, IsEmpty());
  EXPECT_THAT(m_limited, IsEmpty());
}

TEST_F(limit_alleles_test, complex1) {
  start(1);
  add(1, 3, tseq("Block"), 5);
  add(1, 2, tseq("First"), 1);
  add(1, 2, tseq("FirstConflict"), 2);
  add(2, 3, tseq("Second"), 3);
  add(2, 3, tseq("SecondConflict"), 4);
  flush();

  EXPECT_THAT(m_sorts,
              ElementsAre(UnorderedElementsAre(tseq("First"), tseq("FirstConflict"), tseq("Second"),
                                               tseq("SecondConflict"), tseq("Block"))));
  EXPECT_THAT(m_limited,
              UnorderedElementsAre(tseq("FirstConflict"), tseq("SecondConflict"), tseq("Block")));
}
TEST_F(limit_alleles_test, complex2) {
  start(2);
  add(1, 3, tseq("Block"), 5);
  add(1, 2, tseq("First"), 1);
  add(1, 2, tseq("FirstConflict"), 2);
  add(2, 3, tseq("Second"), 3);
  add(2, 3, tseq("SecondConflict"), 4);
  flush();

  EXPECT_THAT(m_sorts,
              ElementsAre(UnorderedElementsAre(tseq("First"), tseq("FirstConflict"), tseq("Second"),
                                               tseq("SecondConflict"), tseq("Block"))));
  EXPECT_THAT(m_limited, UnorderedElementsAre(tseq("Block")));
}

TEST_F(limit_alleles_test, complex3) {
  start(2);
  add(1, 3, tseq("Block"), 3);
  add(1, 2, tseq("First"), 1);
  add(1, 2, tseq("FirstConflict"), 4);
  add(2, 3, tseq("Second"), 2);
  add(2, 3, tseq("SecondConflict"), 5);
  flush();

  EXPECT_THAT(m_sorts,
              ElementsAre(UnorderedElementsAre(tseq("First"), tseq("FirstConflict"), tseq("Second"),
                                               tseq("SecondConflict"), tseq("Block"))));
  EXPECT_THAT(m_limited, UnorderedElementsAre(tseq("FirstConflict"), tseq("SecondConflict")));
}

TEST_F(limit_alleles_test, complex4) {
  start(2);
  add(1, 3, tseq("Block"), 3);
  add(1, 2, tseq("First"), 1);
  add(1, 2, tseq("FirstConflict"), 2);
  add(2, 3, tseq("Second"), 4);
  add(2, 3, tseq("SecondConflict"), 5);
  flush();

  EXPECT_THAT(m_sorts,
              ElementsAre(UnorderedElementsAre(tseq("First"), tseq("FirstConflict"), tseq("Second"),
                                               tseq("SecondConflict"), tseq("Block"))));
  EXPECT_THAT(m_limited, UnorderedElementsAre(tseq("Block")));
}

TEST_F(limit_alleles_test, complex5) {
  start(2);
  add(1, 3, tseq("Block"), 3);
  add(1, 2, tseq("First"), 4);
  add(1, 2, tseq("FirstConflict"), 5);
  add(2, 3, tseq("Second"), 1);
  add(2, 3, tseq("SecondConflict"), 2);
  flush();

  EXPECT_THAT(m_sorts,
              ElementsAre(UnorderedElementsAre(tseq("First"), tseq("FirstConflict"), tseq("Second"),
                                               tseq("SecondConflict"), tseq("Block"))));
  EXPECT_THAT(m_limited, UnorderedElementsAre(tseq("Block")));
}

TEST_F(limit_alleles_test, inserts) {
  start(2);
  add(1, 3, tseq("Block"), 4);
  add(1, 2, tseq("First"), 5);
  add(1, 2, tseq("FirstConflict"), 6);
  add(2, 2, tseq("Insert"), 3);
  add(2, 3, tseq("Second"), 1);
  add(2, 3, tseq("SecondConflict"), 2);
  flush();

  EXPECT_THAT(m_sorts, ElementsAre(UnorderedElementsAre(tseq("First"), tseq("FirstConflict"),
                                                        tseq("Second"), tseq("SecondConflict"),
                                                        tseq("Block"), tseq("Insert"))));
  EXPECT_THAT(m_limited, UnorderedElementsAre(tseq("Block")));
}

TEST_F(limit_alleles_test, inserts2) {
  start(2);
  add(1, 3, tseq("Block"), 4);
  add(1, 2, tseq("First"), 5);
  add(1, 2, tseq("FirstConflict"), 6);
  add(2, 2, tseq("Insert"), 3);
  add(2, 2, tseq("Insert2"), 7);
  add(2, 2, tseq("InsertConflict"), 8);
  add(2, 3, tseq("Second"), 1);
  add(2, 3, tseq("SecondConflict"), 2);
  flush();

  EXPECT_THAT(m_sorts,
              ElementsAre(UnorderedElementsAre(
                  tseq("First"), tseq("FirstConflict"), tseq("Second"), tseq("SecondConflict"),
                  tseq("Block"), tseq("Insert"), tseq("InsertConflict"), tseq("Insert2"))));
  EXPECT_THAT(m_limited, UnorderedElementsAre(tseq("Block"), tseq("InsertConflict")));
}

}  // namespace variants

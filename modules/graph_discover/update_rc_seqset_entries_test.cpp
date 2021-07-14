#include "modules/graph_discover/update_rc_seqset_entries.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class update_rc_seqset_entries_test : public assemble_test {
 public:
  void start() {
      std::unique_ptr<update_rc_seqset_entries> update = make_unique<update_rc_seqset_entries>(m_options, test_output());
      m_update->enable_self_test();
  }
  void flush() {
    m_update->flush();
      EXPECT_TRUE(m_update->self_test_succeeded());
      m_update.reset();
  }
  void run_pass() {
    std::vector<assembly> in_asms = std::move(m_assemblies);
    m_assemblies.clear();

    start();
    for (const auto& a : in_asms) {
      m_update->add(make_unique<assembly>(a));
    }
  }
};

TEST_F(update_rc_seqset_entries_test, simple) {
  use_ref_parts({{100, tseq("abcdefg")}});

  run_update_rc_seqset_entries(0, 1000);

  EXPECT_THAT(m_ref_assemblies, ElementsAre(RefAssemblyIs(100, 100 + tseq("abcdefg").size())));
}

TEST_F(update_rc_seqset_entries_test, part) {
  use_ref_parts({{100, tseq("abcdefg")}});

  run_update_rc_seqset_entries(100 + tseq("a").size(), 100 + tseq("abc").size());

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(100 + tseq("a").size(), 100 + tseq("abc").size())));
}

TEST_F(update_rc_seqset_entries_test, multi_extents) {
  use_ref_parts({{1000, tseq("abc")}, {2000, tseq("def")}, {3000, tseq("ghi")}});

  run_update_rc_seqset_entries(1000 + tseq("a").size(), 3000 + tseq("g").size());
  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(1000 + tseq("a").size(), 1000 + tseq("abc").size()),
                          RefAssemblyIs(2000, 2000 + tseq("def").size()),
                          RefAssemblyIs(3000, 3000 + tseq("g").size())));
}

}  // namespace variants

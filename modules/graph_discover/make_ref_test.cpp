#include "modules/graph_discover/make_ref.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class make_ref_test : public assemble_test {
 public:
  void run_make_ref(aoffset_t start, aoffset_t end, aoffset_t max_chunk_size = 0) {
    {
      auto asms = make_ref_assemblies(m_scaffold, start, end, max_chunk_size);
      auto out = test_output();
      for (auto& a : asms) {
        out->add(std::move(a));
      }
      out->flush();
    }
    EXPECT_THAT(m_non_ref_assemblies, IsEmpty());
  }
};

TEST_F(make_ref_test, simple) {
  use_ref_parts({{100, tseq("abcdefg")}});

  run_make_ref(0, 1000);

  EXPECT_THAT(m_ref_assemblies, ElementsAre(RefAssemblyIs(100, 100 + tseq("abcdefg").size())));
}

TEST_F(make_ref_test, part) {
  use_ref_parts({{100, tseq("abcdefg")}});

  run_make_ref(100 + tseq("a").size(), 100 + tseq("abc").size());

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(100 + tseq("a").size(), 100 + tseq("abc").size())));
}

TEST_F(make_ref_test, multi_extents) {
  use_ref_parts({{1000, tseq("abc")}, {2000, tseq("def")}, {3000, tseq("ghi")}});

  run_make_ref(1000 + tseq("a").size(), 3000 + tseq("g").size());
  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(RefAssemblyIs(1000 + tseq("a").size(), 1000 + tseq("abc").size()),
                          RefAssemblyIs(2000, 2000 + tseq("def").size()),
                          RefAssemblyIs(3000, 3000 + tseq("g").size())));
}

TEST_F(make_ref_test, chunks) {
  use_ref_parts({{100, tseq("abcdef")}});

  run_make_ref(100 + tseq("a").size(), 100 + tseq("abc").size(), tseq("a").size());

  EXPECT_THAT(m_ref_assemblies,
              ElementsAre(
                          RefAssemblyIs(100 + tseq("a").size(), 100 + tseq("ab").size()),
                          RefAssemblyIs(100 + tseq("ab").size(), 100 + tseq("abc").size())));
}

}  // namespace variants

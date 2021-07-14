#include "modules/variants/anchor_split.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class anchor_split_test : public assemble_test {
 public:
  void SetUp() { m_options.trace_reference_assemblies = true;
    m_anchor_split.emplace(m_options, test_output()); }
  void add(assembly a) { m_anchor_split->add(make_unique<assembly>(a)); }
  void flush() { m_anchor_split.reset(); }

 protected:
  boost::optional<anchor_splitter> m_anchor_split;
};

TEST_F(anchor_split_test, splits) {
  assembly a;
  a.seq = tseq("abcdef");
  a.left_offset = 100;
  a.right_offset = 200;
  a.left_anchor_len = tseq("ab").size();
  a.right_anchor_len = tseq("f").size();

  add(std::move(a));
  flush();

  EXPECT_THAT(m_assemblies, UnorderedElementsAre(
                                RefAssemblyIs(100, 100 + tseq("ab").size()),
                                AssemblyIs(100 + tseq("ab").size(), tseq("cde"),
                                           200 - tseq("f").size()),
                                RefAssemblyIs(200 - tseq("f").size(), 200)));
}

}  // namespace variants

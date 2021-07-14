#include "modules/graph_discover/graph_trim_ref.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/assemble_testutil.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class graph_trim_ref_test : public assemble_test {
 public:
  graph_trim_ref_test() { m_graph_trim_ref.emplace(m_options, test_output()); }
  void add(assembly a) { m_graph_trim_ref->add(make_unique<assembly>(a)); }
  void flush() {
    m_graph_trim_ref->flush();
    m_graph_trim_ref.reset();
  }

  void add_ref(aoffset_t left_offset, dna_sequence seq) {
    assembly a;
    a.matches_reference = true;
    a.left_offset = left_offset;
    a.right_offset = left_offset + seq.size();
    a.seq = seq;

    add(a);
  }

  void add_var(optional_aoffset left_offset, dna_sequence seq, optional_aoffset right_offset) {
    assembly a;
    a.left_offset = left_offset;
    a.seq = seq;
    a.right_offset = right_offset;
    add(a);
  }

 protected:
  boost::optional<graph_trim_ref> m_graph_trim_ref;
};

TEST_F(graph_trim_ref_test, simple_ref) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  add_ref(0, tseq("abcd"));
  flush();

  EXPECT_THAT(m_assemblies, ElementsAre(RefAssemblyIs(0, tseq("abcd").size())));
}

TEST_F(graph_trim_ref_test, trim_left_anchor) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  add_ref(0, tseq("abcd"));
  add_var(0, tseq("ab") + dna_G, optional_aoffset::none);
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(RefAssemblyIs(0, tseq("ab").size()),
                                   RefAssemblyIs(tseq("ab").size(), tseq("abcd").size()),
                                   AssemblyIs(tseq("ab").size(), dna_G, optional_aoffset::none)));
}

TEST_F(graph_trim_ref_test, trim_right_anchor) {
  use_ref_parts({{0, tseq("abcdefghijklmnopqrstuvwxyz")}});
  add_ref(0, tseq("abcd"));
  add_var(optional_aoffset::none, dna_G + tseq("cd"), tseq("abcd").size());
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(RefAssemblyIs(0, tseq("ab").size()),
                                   RefAssemblyIs(tseq("ab").size(), tseq("abcd").size()),
                                   AssemblyIs(optional_aoffset::none, dna_G, tseq("ab").size())));
}

TEST_F(graph_trim_ref_test, trim_both) {
  use_ref_parts({{0, tseq("abcdef")}});
  add_ref(0, tseq("abcdef"));
  add_var(tseq("ab").size(), tseq("c") + dna_G + tseq("e"), tseq("abcde").size());
  flush();

  EXPECT_THAT(m_assemblies,
              UnorderedElementsAre(RefAssemblyIs(0, tseq("abc").size()),
                                   RefAssemblyIs(tseq("abc").size(), tseq("abcd").size()),
                                   RefAssemblyIs(tseq("abcd").size(), tseq("abcdef").size()),
                                   AssemblyIs(tseq("abc").size(), dna_G, tseq("abcd").size())));
}

TEST_F(graph_trim_ref_test, trim_to_ref) {
  use_ref_parts({{0, tseq("abcdef")}});
  add_ref(0, tseq("abcdef"));
  add_var(tseq("abc").size(), tseq("d"), tseq("abcd").size());
  flush();

  EXPECT_THAT(m_assemblies, UnorderedElementsAre(RefAssemblyIs(0, tseq("abcdef").size())));
}

}  // namespace variants

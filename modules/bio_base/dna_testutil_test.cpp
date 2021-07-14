#include "modules/bio_base/dna_testutil.h"

#include <gtest/gtest.h>

TEST(dna_testutils, enable_disable_printer) {
  using namespace dna_testutil;

  dna_sequence seq = tseq("hi!");
  dna_slice slice = seq;

  // By default, dna_testutil_env should enable expansion.
  {
    std::stringstream os;
    std::string expected_hi =
        "\"" + dna_test_sequence("hi!").as_string() + "\" (tseq(\"hi!\"))";

    os << seq;
    os << slice;

    EXPECT_EQ(expected_hi + expected_hi, os.str());
  }

  // But we should be able to disable it.
  disable_test_sequence_expansion();

  {
    std::stringstream os;

    dna_sequence seq = tseq("hi!");
    dna_slice slice = seq;

    os << seq;
    os << slice;

    EXPECT_EQ(tseq("hi!hi!"), os.str());
  }

  // And re-enable it.
  enable_test_sequence_expansion();

  {
    std::stringstream os;
    std::string expected_hi =
        "\"" + dna_test_sequence("hi!").as_string() + "\" (tseq(\"hi!\"))";

    os << seq;
    os << slice;

    EXPECT_EQ(expected_hi + expected_hi, os.str());
  }
}

TEST(dna_testutils, dna_test_sequence) {
  EXPECT_EQ(dna_sequence("CAAAAAAAAC"), dna_test_sequence(std::string({0})));
  EXPECT_EQ(dna_sequence("CTAAAAAAAC"), dna_test_sequence(std::string({1})));
  EXPECT_EQ(dna_sequence("CTTTTTTTAC"), dna_test_sequence(std::string({127})));
  EXPECT_EQ(dna_test_sequence("hello there"),
            dna_test_sequence("h") + dna_test_sequence("ello there"));
  EXPECT_EQ(k_dna_test_sequence_length, dna_test_sequence("x").size());
}

TEST(dna_testutils, debug_print) {
  EXPECT_EQ("\n  \"" + dna_test_sequence("hi!").as_string() +
                "\" (tseq(\"hi!\"))",
            ::testing::PrintToString(dna_test_sequence("hi!")));
  EXPECT_EQ("\n  \"" + dna_test_sequence("hi!").rev_comp().as_string() +
                "\" (tseq_rc(\"hi!\"))",
            ::testing::PrintToString(dna_test_sequence("hi!").rev_comp()));
  EXPECT_EQ("\n  \"GATTACA" + dna_test_sequence("x").as_string() +
                "TAGGA\" (\"GATTACA\" + tseq(\"x\") + \"TAGGA\")",
            ::testing::PrintToString(dna_sequence("GATTACA") +
                                     dna_test_sequence("x") +
                                     dna_sequence("TAGGA")));
  EXPECT_EQ("\n  \"\"", ::testing::PrintToString(dna_test_sequence("")));
  EXPECT_EQ("\n  \"CAAACAAAA" + dna_test_sequence("q").as_string() +
                "\" (\"CAAACAAAA\" + tseq(\"q\"))",
            ::testing::PrintToString(dna_sequence("CAAACAAAA") +
                                     dna_test_sequence("q")));
}

TEST(dna_testutils, short_aliases) {
  using dna_testutil::tseq;
  using dna_testutil::tseq_rc;
  EXPECT_EQ(dna_test_sequence("howdy"), tseq("howdy"));
  EXPECT_EQ(dna_test_sequence("howdy").rev_comp(), tseq_rc("howdy"));
}

TEST(dna_testutils, constants) {
  using namespace dna_testutil;
#define CHECK_BASE(base) EXPECT_EQ(dna_sequence(#base), dna_##base)
  CHECK_BASE(A);
  CHECK_BASE(C);
  CHECK_BASE(G);
  CHECK_BASE(T);
#undef CHECK_BASE
}

TEST(dna_testutils, drop) {
  using dna_testutil::drop_front;
  using dna_testutil::drop_back;

  EXPECT_EQ(drop_front(2, dna_sequence("ACTGA")), dna_sequence("TGA"));
  EXPECT_EQ(drop_front(5, dna_sequence("ACTGA")), dna_sequence(""));
  EXPECT_EQ(drop_back(2, dna_sequence("ACTGA")), dna_sequence("ACT"));
  EXPECT_EQ(drop_back(5, dna_sequence("ACTGA")), dna_sequence(""));
}

#include "modules/bio_base/sv_call_testutil.h"
#include "modules/bio_base/dna_testutil.h"

using namespace testing;
using namespace dna_testutil;
using namespace sv_call_testutil;

TEST(sv_call_testutil_test, sv_call_matcher) {
  sv_call call;
  call.position.scaffold_id = 3;
  call.position.position = 73;

  allele a;
  a.seq = tseq("a");
  call.alleles.push_back(a);

  allele b;
  b.seq = tseq("b");
  call.alleles.push_back(b);

  EXPECT_THAT(call, SvCall(SeqPosition(3, 73),
                           ElementsAre(
                               Allele(tseq("a")),
                               Allele(tseq("b")))));
  EXPECT_THAT(call, Not(SvCall(SeqPosition(3, 74),
                           ElementsAre(
                               Allele(tseq("a")),
                               Allele(tseq("b"))))));
  EXPECT_THAT(call, Not(SvCall(SeqPosition(4, 73),
                           ElementsAre(
                               Allele(tseq("a")),
                               Allele(tseq("b"))))));
  EXPECT_THAT(call, Not(SvCall(SeqPosition(3, 73),
                           ElementsAre(
                               Allele(tseq("a"))))));
  EXPECT_THAT(call, Not(SvCall(SeqPosition(3, 73),
                           ElementsAre(
                               Allele(tseq("a")),
                               Allele(tseq("c"))))));
}

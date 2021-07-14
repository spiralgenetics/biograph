#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>

#include "modules/variants/assemble_testutil.h"
#include "modules/variants/big_assemble_testutil.h"

using namespace testing;
using namespace variants;

constexpr bool k_full_human = false;

class tinyhuman_test : public big_assemble_test {
 public:
  void SetUp() override {
    if (k_full_human) {
      use_biograph("HG001_frc_8.bg");
    } else {
      use_biograph("tinyhuman-loop.bg");
    }
    // TODO(nils): Figure out why we need so many:
    m_options.max_coverage_paths = 20;
    m_options.use_bidir_tracer = true;
  }
};

// 2       220863241       rs12162365      A       T       50      PASS
// platforms=4;platformnames=Illumina,CG,10X,Solid;datasets=5;
// datasetnames=HiSeqPE300x,CGnormal,10XChromium,SolidPE50x50bp,SolidSE75bp;
// callsets=6;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,
// 10XGATKhaplo,SolidPE50x50GATKHC,SolidSE75GATKHC;datasetsmissingcall=IonExome;
// callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,
// CS_HiSeqPE300xfreebayes_callable
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|0:614:152,134:167,155:99:0/1:.:PATMAT
TEST_F(tinyhuman_test, rs12162365) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("2", "220863241", "A", "T", "0/1");
  // call_at("2", "220863241", 200, 200);

  // EXPECT_EQ(dna_base('A'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "A")));         // reference
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "T")));  // variant
}

// 13 44617398        rs4562968       A       C       50      PASS
// platforms=4;platformnames=Illumina,CG,10X,Solid;datasets=5;datasetnames=HiSeqPE300x,CGnormal,10XChromium,SolidPE50x50bp,SolidSE75bp;callsets=6;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo,SolidPE50x50GATKHC,SolidSE75GATKHC;datasetsmissingcall=IonExome;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable,CS_10XGATKhaplo_callable;filt=CS_SolidPE50x50GATKHC_filt,CS_SolidSE75GATKHC_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|1:493:0,196:57,241:99:1/1:.:PATMAT
TEST_F(tinyhuman_test, rs4562968) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("13", "44617398", "A", "C", "1/1");
  // call_at("13", "44617398", 400, 400);
  // // 70492481 from reverse

  // EXPECT_EQ(dna_base('A'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // // homozygous variant; no reference allele.
  // EXPECT_THAT(m_assemblies, Not(Contains(RefAt(m_call_pos, "A"))));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "C")));
}

// 1       187384691       rs7527494       A       G       50      PASS
// platforms=4;platformnames=Illumina,CG,10X,Solid;datasets=5;datasetnames=HiSeqPE300x,CGnormal,10XChromium,SolidPE50x50bp,SolidSE75bp;callsets=6;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo,SolidPE50x50GATKHC,SolidSE75GATKHC;datasetsmissingcall=IonExome;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable;filt=CS_HiSeqPE300xfreebayes_filt,CS_SolidSE75GATKHC_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|0:626:159,138:168,161:99:0/1:.:PATMAT
TEST_F(tinyhuman_test, rs7527494) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("1", "187384691", "A", "G", "0/1");
  // call_at("1", "187384691", 400, 350);

  // // Assembly id=59 [187384690:0, 187384691:0), right pair count=0, left pair count = 2, !ref: G
  // // (len=1)

  // EXPECT_EQ(dna_base('A'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);

  // EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "A")));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "G")));
}

// 8  15978063        rs4338104       C       T       50      PASS
// platforms=5;platformnames=Illumina,CG,10X,Ion,Solid;datasets=5;datasetnames=HiSeqPE300x,CGnormal,10XChromium,IonExome,SolidSE75bp;callsets=6;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo,IonExomeTVC,SolidSE75GATKHC;datasetsmissingcall=SolidPE50x50bp;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable,CS_10XGATKhaplo_callable,CS_IonExomeTVC_callable;filt=CS_HiSeqPE300xGATK_filt,CS_HiSeqPE300xfreebayes_filt,CS_SolidSE75GATKHC_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|1:551:0,124:36,36:99:1/1:.:PATMAT
TEST_F(tinyhuman_test, rs4338104) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("8", "15978063", "C", "T", "1/1");
  // call_at("8", "15978063", 250, 250);

  // // Not enough overlap reads for the section of repeated As from
  // // 8:15978126 to 8:15978150.

  // EXPECT_EQ(dna_base('C'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // // homozygous variant; should not contain reference.
  // EXPECT_THAT(m_assemblies, Not(Contains(RefAt(m_call_pos, "A"))));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "T")));
}

// 1       169650747       rs7543304       T       C       50      PASS
// platforms=4;platformnames=Illumina,CG,10X,Solid;datasets=4;datasetnames=HiSeqPE300x,CGnormal,10XChromium,SolidSE75bp;callsets=5;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo,SolidSE75GATKHC;datasetsmissingcall=IonExome,SolidPE50x50bp;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable,CS_10XGATKhaplo_callable
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|0:718:149,168:197,207:99:0/1:.:PATMAT
TEST_F(tinyhuman_test, rs7543304) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("1", "169650747", "T", "C", "0/1");
  // call_at("1", "169650747", 300, 300);

  // EXPECT_EQ(dna_base('T'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // if (k_test_genotyping) {
  //   EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "T")));
  // }
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "C")));
}

// 8       94243253        rs57983550      AC      A       50      PASS
// platforms=3;platformnames=Illumina,CG,10X;datasets=3;datasetnames=HiSeqPE300x,CGnormal,10XChromium;callsets=4;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo;datasetsmissingcall=IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable;filt=CS_HiSeqPE300xGATK_filt,CS_HiSeqPE300xfreebayes_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS 0|1:396:97,50:28,36:99:0|1:94243231_C_G:PATMAT
TEST_F(tinyhuman_test, rs57983550) {
  SCOPED_BIG_ASM_TEST();

  // On the left, we need to start at 94243169, which is a branch off at A.
  // Then, at 94243231, we see a "G" on one allele and the other stays on
  // reference at "C".
  // The allele that branched off to "G" should see a delete at 94243253.
  // It should notice this when it rejoins the read from 94243273 to 94243365.
  //
  //
  // As of 2018-03-14, we were getting:
  //  , sv_call at scaffold 7 position 94243304 with 3 alleles:
  // allele 0: TTGGGTGATAGAATGAGACTCTGTCTCAAAAAAAAA depth: {5, 5, 5, 5, 5, 5, 5,
  // 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 2, 2, 1, 1, 1,
  // 1, 1, 2, 2, 2, 2}
  // allele 1: CATTGCACTCCAGCCTGGGTGACAGAGAGAGACTCTGTC ids: {2*2} depth: {1, 1,
  // 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  // 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
  // allele 2: TTGGGTGATAGAATGAGAAACCCTTTCCAAAAAAAA ids: {4*2} depth: {1, 1, 1,
  // 1, 1, 1, 1, 1, 1, 1, 1, 1}
  //
  // allele 1 needs to get filtered.  but other than that, it's
  // getting called as a change from "A" to "C" instead of an insert
  // of the "C".

  /*  m_trace_enabled = true;

  add_trace(
      94243020 + 148, 94243342 - 88,
      "AAAAATACAAAAAATTAGCTGGGTGTGGTGGCAGGCACCTGTAATCCCAGCTACTCAGGAGGGTGAGGCATGAGAATCATTTGAA"
      // rc: TTCAAATGATTCTCATGCCTCACCCTCCTGAGTAGCTGGGATTACAGGTGCCTGCCACCACACCCAGCTAATTTTTTGTATTTTT

      ); 
  */

  run_vcf_test("8", "94243253", "AC", "A", "0/1");
  // call_at("8", "94243254", 100, 100);

  // EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "C")));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "")));
}

// 16      24600517        rs7200722       T       C       50      PASS
// platforms=3;platformnames=Illumina,CG,10X;datasets=3;datasetnames=HiSeqPE300x,CGnormal,10XChromium;callsets=4;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo;datasetsmissingcall=IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|0:578:111,99:172,137:99:0/1:.:PATMAT
TEST_F(tinyhuman_test, rs7200722) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("16", "24600517", "T", "C", "0/1");
  // call_at("16", "24600517", 400, 400);

  // EXPECT_EQ(dna_base('T'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "T")));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "C")));
}

// 14  38766393        rs71433950      C       CTA     50      PASS
// platforms=1;platformnames=Illumina;datasets=1;datasetnames=HiSeqPE300x;callsets=2;callsetnames=HiSeqPE300xGATK,HiSeqPE300xfreebayes;datasetsmissingcall=CGnormal,10XChromium,IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable;difficultregion=AllRepeats_lt51bp_gt95identity_merged_slop5
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    0/1:514:105,124:105,124:99:0/1
// end of scaffold is 107289540
// 1072895340 - 38766392 = 68523148

TEST_F(tinyhuman_test, rs71433950) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("14", "38766393", "C", "CTA", "0/1");
  // call_at("14", "38766393", 3, 3);

  // auto vcf_added = m_call_ref_it;
  // EXPECT_EQ(dna_base('C'), *vcf_added) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                      << dna_slice(m_call_ref_it, 10);
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos + 1, 0, "TA")));
}

// 4     14142467        rs7694979       G       T       50      PASS
// platforms=3;platformnames=Illumina,10X,CG;datasets=3;datasetnames=HiSeqPE300x,10XChromium,CGnormal;callsets=4;callsetnames=HiSeqPE300xGATK,HiSeqPE300xfreebayes,10XGATKhaplo,CGnormal;datasetsmissingcall=IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable,CS_HiSeqPE300xfreebayes_callable,CS_10XGATKhaplo_callable;filt=CS_HiSeqPE300xGATK_filt,CS_SolidSE75GATKHC_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|0:651:179,167:0,0:99:0/1:.:PATMAT
TEST_F(tinyhuman_test, rs7694979) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("4", "14142467", "G", "T", "0/1");
  // call_at("4", "14142467", 200, 200);

  // // 191044276 14142467 - p
  // // 176901809

  // // Not enough overlap reads for the section of repeated As from
  // // 8:15978126 to 8:15978150.

  // EXPECT_EQ(dna_base('G'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "G")));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "T")));
}

// 19      33527793        rs11879039      T       C       50      PASS
// platforms=3;platformnames=Illumina,CG,10X;datasets=3;datasetnames=HiSeqPE300x,CGnormal,10XChromium;callsets=4;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo;datasetsmissingcall=IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable;filt=CS_10XGATKhaplo_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS 1|0:469:125,88:146,99:99:0|1:33527784_G_C:PATMAT
TEST_F(tinyhuman_test, rs11879039) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("19", "33527793", "T", "C", "0/1");
  // call_at("19", "33527793", 1700, 1700);

  // EXPECT_EQ(dna_base('T'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "T")));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "C")));
}

// 7 69271250 rs7809145 G A 50 PASS
// platforms=4;platformnames=Illumina,CG,10X,Solid;datasets=4;datasetnames=HiSeqPE300x,CGnormal,10XChromium,SolidSE75bp;callsets=5;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo,SolidSE75GATKHC;datasetsmissingcall=IonExome,SolidPE50x50bp;callable=CS_HiSeqPE300xGATK_callable,CS_CGnormal_callable,CS_HiSeqPE300xfreebayes_callable,CS_10XGATKhaplo_callable;filt=CS_HiSeqPE300xfreebayes_filt,CS_10XGATKhaplo_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS 1|0:452:68,92:106,122:99:0/1:.:PATMAT
TEST_F(tinyhuman_test, rs7809145) {
  // 1241-12512 forward
  // Rev: 159128663 - 69271249 = 89857414
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("7", "69271250", "G", "A", "0/1");
  // call_at("7", "69271250", 400, 400);

  // EXPECT_EQ(dna_base('G'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // if (k_test_genotyping) {
  //   EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "G")));
  // }
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "A")));
}

#if 0
// False positive: 13      108575772       .       T       G       100     PASS
// NS=1;DP=792;AID=35483,35784,35728,35780,35687,35643,35435,35594,35785;FW=0.5;BQ=30
// GT:DP:AD:OV     0/1:792:82,709:70
TEST_F(tinyhuman_test, fp_13_108575772) {
  call_at("13", "108575772", 200, 200);

  EXPECT_EQ(dna_base('T'), *m_call_ref_it);
  EXPECT_THAT(
      m_assemblies,
      Not(Contains(SvCall(m_call_pos, ElementsAre(Allele("T"),  // reference
                                                  Allele("G")   // variant
                                                  )))));
}
#endif

// 16      58817872        rs147377264     A       C       50      PASS
// platforms=2;platformnames=Illumina,10X;datasets=2;datasetnames=HiSeqPE300x,10XChromium;callsets=3;callsetnames=HiSeqPE300xGATK,HiSeqPE300xfreebayes,10XGATKhaplo;datasetsmissingcall=CGnormal,IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable,CS_HiSeqPE300xfreebayes_callable
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS    1|1:544:0,277:0,277:99:1/1:.:PATMAT

TEST_F(tinyhuman_test, rs147377264) {
  SCOPED_BIG_ASM_TEST();

  /*
    m_trace_enabled = true;
    add_trace(58817704 + 148, 58818036 - 148, "GATGTAGTTTCTTCCTAGCCTCGATGGTCTTTACGT");
    add_trace(58817704 + 148, 58817952 - 65,  "GATGTAGTTTCTTCCTAGCCTCGATGGTCTTTACG");
  */

  run_vcf_test("16", "58817872", "A", "C", "1/1");
// call_at("16", "58817872", 200, 200);
// EXPECT_EQ(dna_base('A'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
//                                          << dna_slice(m_call_ref_it, 10);
// // homozygous variant should not have reference present.
// // TODO(nils): Uncomment when we have genotyping integrated:
// // EXPECT_THAT(m_assemblies, Not(Contains(RefAt(m_call_pos, "A"))));
// EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "C")));
#if 0
  // Don't want these:
  // 16	58817871	.	CATCG	TCTTC	100	PASS
  // NS=1;DP=118;AID=197024,190193,197172;FW=0.4;BQ=30	GT:DP:AD:OV
  // ./.:118:0,54,60,3:70
  // 16	58817875	.	G	A	100	PASS
  // NS=1;DP=118;AID=197024,190193,197172;FW=0.49;BQ=30	GT:DP:AD:OV
  // ./.:118:0,54,60,3:70
  // 16	58817882	.	T	G	100	PASS	NS=1;DP=54;AID=197172;FW=0.5;BQ=30
  // GT:DP:AD:OV	1/1:54:0,54:70
  EXPECT_THAT(m_assemblies,
              Not(Contains(SvCall(m_ref->get_seq_position(m_flat_call_pos - 1),
                                  ElementsAre(Allele("CATCG"),  // reference
                                              Allele("TCTTC")   // variannt
                                              )))));
  EXPECT_THAT(m_assemblies,
              Not(Contains(SvCall(m_ref->get_seq_position(m_flat_call_pos + 3),
                                  ElementsAre(Allele("G"),  // reference
                                              Allele("A")   // variannt
                                              )))));
  EXPECT_THAT(m_assemblies,
              Not(Contains(SvCall(m_ref->get_seq_position(m_flat_call_pos + 10),
                                  ElementsAre(Allele("T"),  // reference
                                              Allele("G")   // variannt
                                              )))));
#endif
}

// 17 6974488 rs71383454 T TC 50 PASS
// platforms=3;platformnames=Illumina,CG,10X;datasets=3;datasetnames=HiSeqPE300x,CGnormal,10XChromium;callsets=4;callsetnames=HiSeqPE300xGATK,CGnormal,HiSeqPE300xfreebayes,10XGATKhaplo;datasetsmissingcall=IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable;filt=CS_CGnormal_filt;difficultregion=AllRepeats_lt51bp_gt95identity_merged_slop5,SimpleRepeat_imperfecthomopolgt10_slop5
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS 0|1:443:104,144:68,103:99:0/1:.:PATMAT
TEST_F(tinyhuman_test, rs71383454) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("17", "6974488", "T", "TC", "0/1");
  // call_at("17", "6974488", 10, 10);

  // EXPECT_EQ(dna_base('T'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos + 1, 0, "C")));
}

// 1 38362435 rs182445103 C G 50 PASS
// platforms=1;platformnames=Illumina;datasets=1;datasetnames=HiSeqPE300x;callsets=2;callsetnames=HiSeqPE300xGATK,HiSeqPE300xfreebayes;datasetsmissingcall=CGnormal,10XChromium,IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable;filt=CS_HiSeqPE300xfreebayes_filt;difficultregion=AllRepeats_lt51bp_gt95identity_merged_slop5
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS
// 0|1:448:104,97:104,97:99:1|0:38362377_G_A:PATMAT
TEST_F(tinyhuman_test, bidir_maybe_rs182445103) {
  if (!m_options.use_bidir_tracer) {
    return;
  }
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("1", "38362435", "C", "G", "0/1");
  // call_at("1", "38362435", 50, 50);
  // // Scaffold end pos: 249240621
  // // Call distance from end: 210878187

  // EXPECT_EQ(dna_base('C'), *m_call_ref_it) << dna_slice(m_call_ref_it - 10, 10) << " "
  //                                          << dna_slice(m_call_ref_it, 10);

  // EXPECT_THAT(m_assemblies, Contains(RefAt(m_call_pos, "C")));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos, 1, "G")));
}

// Make sure this assembly aligns correctly:
// 1 7163493 rs151171782 C CACA 50 PASS
// platforms=2;platformnames=Illumina,10X;datasets=2;datasetnames=HiSeqPE300x,10XChromium;callsets=3;callsetnames=HiSeqPE300xGATK,HiSeqPE300xfreebayes,10XGATKhaplo;datasetsmissingcall=CGnormal,IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable,CS_HiSeqPE300xfreebayes_callable;filt=CS_10XGATKhaplo_filt
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS 1|1:489:0,238:0,238:198:1/1:.:PATMAT
//
// It should not align as:
// 1 7163495 .  C CAAC 100 PASS NS=1;AID=40619,54054
// PI:GT:OV:PG:DP:AD:GQ 40619:0/1:124:0|1:39:0,39:99

// CACT
// CACAACT
// CACAACT

// The change is from:
// ...CACT...
// to:
// ...CACAACT
//
// Normalization says this should be "C" + "ACA" + "AACT",
// not "CAC" + "AAC" + "T".

TEST_F(tinyhuman_test, rs151171782) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("1", "7163493", "C", "CACA", "1/1");
  // call_at("1", "7163493", 50, 50);
  // // Scaffold end pos: 249240621
  // // Call distance from end: 210878187

  // EXPECT_EQ(dna_base('C'), *m_call_ref_it) << dna_slice(m_call_ref_it - 100, 100) << " "
  //                                          << dna_slice(m_call_ref_it, 100);

  // EXPECT_THAT(m_assemblies, Not(Contains(RefAt(m_call_pos, "C"))));
  // EXPECT_THAT(m_assemblies, Contains(VariantAt(m_call_pos + 1, 0, "ACA")));
}

// False positive:
// 7       87583817        .       T       C       100     PASS    NS=1;AID=11453891,11453928
// PI:GT:OV:PG:DP:AD:GQ    11453891:0/1:113:0|1:634:235,399:62.934
//
// This occurs in a region that's similar to a lot of other places in
// the human genome, so we need to make sure we're tracing the path
// for the right place.
TEST_F(tinyhuman_test, DISABLED_fp_7_87583817) {
  call_at("7", "87583817", 200, 200);

  EXPECT_EQ(dna_base('T'), *m_call_ref_it) << dna_slice(m_call_ref_it - 100, 100) << " "
                                           << dna_slice(m_call_ref_it, 100);

  EXPECT_THAT(m_assemblies, Not(Contains(VariantAt(m_call_pos, 1, "C"))));
}

// 1 47634413 .  G C 100 PASS NS=1;AID=388745,388730,388736,388729
// PI:GT:OV:PG:DP:AD:GQ 388745:0/1:117:0|1:40:0,40:99
TEST_F(tinyhuman_test, fp_1_47634413) {
  m_options.pop_trace_anchor_drop = false;
  call_at("1", "47634413", 50, 50);

  EXPECT_EQ(dna_base('G'), *m_call_ref_it) << dna_slice(m_call_ref_it - 100, 100) << " "
                                           << dna_slice(m_call_ref_it, 100);

  EXPECT_THAT(m_assemblies, Not(Contains(VariantAt(m_call_pos, 1, "C"))));
}

// 11 47943677 .  A C 100 PASS
// NS=1;AID=16285047,16283352,16283389,16285070,16285121,16285112
// PI:GT:OV:PG:DP:AD:GQ:PDP:PAD 16285047:0/1:115:0|1:31:22,9:2:0,2:100
//
// All the evidence for the "C" is on one strand.
TEST_F(tinyhuman_test, DISABLED_fp_11_47943677) {
  call_at("11", "47943677", 50, 50);

  EXPECT_EQ(dna_base('A'), *m_call_ref_it) << dna_slice(m_call_ref_it - 100, 100) << " "
                                           << dna_slice(m_call_ref_it, 100);

  EXPECT_THAT(m_assemblies, Not(Contains(VariantAt(m_call_pos, 1, "C"))));
}

// 14 29194382 .  CCACACA CCA,C 50 PASS
// platforms=1;platformnames=Illumina;datasets=1;datasetnames=HiSeqPE300x;callsets=2;callsetnames=HiSeqPE300xGATK,HiSeqPE300xfreebayes;datasetsmissingcall=CGnormal,10XChromium,IonExome,SolidPE50x50bp,SolidSE75bp;callable=CS_HiSeqPE300xGATK_callable;filt=CS_HiSeqPE300xfreebayes_filt;difficultregion=AllRepeats_lt51bp_gt95identity_merged_slop5
// GT:DP:ADALL:AD:GQ:IGT:IPS:PS
// 1|2:219:0,57,51:0,57,51:99:2/1:.:PATMAT
//
// Unjoined, we're seeing:
// 14      29194382        .       CCACA   C       100     PASS    NS=1;AID=19317290       PI:GT:OV:PG:GQ:DP:AD:PDP:PAD    19317290:0/1:126:0|1:31:16:5,11:1:0,1
// 14      29194382        .       CCACACA C       100     PASS    NS=1;AID=19318954,19317282      PI:GT:OV:PG:GQ:DP:AD:PDP:PAD    19318954:1/1:111:1|1:50:32:5,16:11:0,10
//
// The second genotype should be 0/1, not 1/1.
// TODO(nils): These two should probably be joined into a single test.
TEST_F(tinyhuman_test, fn_14_29194382_1) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("14", "29194382", "CCACA", "C", "0/1");
}

TEST_F(tinyhuman_test, fn_14_29194382_2) {
  SCOPED_BIG_ASM_TEST();

  run_vcf_test("14", "29194382", "CCACACA", "C", "0/1");
}

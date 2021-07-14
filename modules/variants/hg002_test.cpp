#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>

#include "modules/variants/assemble_testutil.h"
#include "modules/variants/big_assemble_testutil.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/pop_tracer.h"

using namespace testing;
using namespace variants;
using namespace dna_testutil;

class hg002_test : public big_assemble_test {
 public:
  void SetUp() override { use_biograph("HG002-NA24385-50x.bg"); }
};

// There's a big deletion here that's supposedly not cought by bagpipe (vcf, so
// 1-based index):
// 4 34779885        sv_680720      ATGGT...GCCCA G
// It looks like the deletion is actually at  [34779892,34828950)?
TEST_F(hg002_test, sv_4_34779885) {
  // At 4:34779884, there's a structural change of 49058 bases.
  // Most of this is a deletion, but there's another change at the end
  // This should end at 34779883 + 49058 = 34828941,
  // plus anchor of +148 = about 34829089
  //  call_region("4", 34779800, 34830000);
  //

  run_vcf_test("4", "34779885", LongSequenceMatches("ATGGT.*GCCCAGGCG", 49059), "A", "1/1");

  // We originally detected the deletion at 34779892, but normalized it's at 34779885.
  // With the vcf padded base, which is A, it starts at 34779884.
  //  call_region("4", 34770000, 34830000);
  //  EXPECT_THAT(m_assemblies, Contains(AllOf(VariantAt(34779884, 49059, "A"), GenotypeIs("1/1"))))
  //      << dump_sv_assemblies(40000);
}

// 2       34695830        HG2_PB_SVrefine2PBcRplusDovetail_5578
// TGTCTTAGCCCAAAATCTCCTTAAGATGATAAGCAA...GAAAAAATACAGCCTGCTGATGTCCTTCAGAAA
// GTAGGGGT        20      PASS
// END=34736567;SVTYPE=DEL;SVLEN=-40730;ClusterIDs=HG3_Ill_MetaSV_2617:HG3_PB_pbsv_1791:Hn
// TODO(nils): Fix this test
TEST_F(hg002_test, DISABLED_sv_2_34695830) {
  run_vcf_test("2", "34695830", LongSequenceMatches("TGTCTTAGCCCAAAATCTCCTTAAGATGATAAGCAA.*GAAAAAATACAGCCTGCTGATGTCCTTCAGAAA", 40738), "GTAGGGGT", "0/1");

  //  call_region("2", 34695678, 34736715);
  //  EXPECT_THAT(m_assemblies,
  //              Contains(AllOf(VariantAt(34695832, 40735, "AGGGGT"), GenotypeIs("0/1"))))
  //      << dump_sv_assemblies(30000);
}

// ref len: 20054
// alt len: 3
// 4:64694329 HG2_PB_SVrefine2Falcon2Bionano_2184
// TATTGTGTATTATCACTAGTT..TATTATGTTATATTATATTATATAAGATATTATATTATATTTATGTATTAATATATTAATATGTTTATAT
// AAA
// end = 64694329 + 20054 = 64714383
//
// Detected raw assembly: [64694180:148+20301=64714481:90), !ref:
// TTTGTCTTCTATTCTTTGCTCTTCTGAAATAGGGTTTTTATTGCTGTTATTTTTTGCCCATGCTTCTTTATTAATAGTATTTCATTATGACAACTCATATTTTAAGTCCATAGATCTATGACCAAGAGAACCACAATTTGAAATCATG
// AAAA TTATAAA G
// ATAATATATAAACATAATATATTGTAAACATAATATATGGTTTATATATAATATATATTATATAATATATATTATATTATAATCATATAT
// id=24 score=3400
//
// At beginning of deletion:
//
// Calling around 64694329
// Ref before: ATTTTAAGTCCATAGATCTATGACCAAGAGAACCACAATTTGAAATCATG
//     after:
//     TATTGTGTATTATCACTAGTTAAAGACAATTAGAAAAGCAAACAAATGTGATTCATGTACAGATCCTGAACCACAGATTCAAGGCTCAATGTCCAAAATGGACAGGACTTTGAGTTGTCTCTCTGAAAATATTATATGTGTATTTTGCATCAGAAAATAATAATAAATTAAATATTTGGTGACTAGAGTGGTAGTTTCCAATAGTCATTGTAACTATTAAAAATGTGATTCTATTCTTTTCCCAGAAAAAGTAAAATTTTACATTTTGATCAAATTTAAACTTGGACATACTAGCATAAA
// Scaffold end pos: 191044276
// Call distance from end: 126349948
//
// At end of deletion:
// Calling around 64714383
// Ref before: TATAAGATATTATATTATATTTATGTATTAATATATTAATATGTTTATAT
//     after:  A TTATAAA C ATAATATATAAACATAATATATTGTAAACATAATATATGGT
// Scaffold end pos: 191044276
// Call distance from end: 126329894
//
// This traces to a dead end that doesn't contain a read that entirely
// matches reference.  However, the dead end has enough overlap with
// reference that it should detect the 20054 bases replaced with
// "AAA".
TEST_F(hg002_test, sv_4_64694329) {
  m_options.min_overlap = 70;
  // add_assembly_trace(246);
  // add_assembly_trace(247);

  // We are aligning this as a SNP of A->T, followed by a 20054 length
  // SV instead of a 20056 length SV.
  //
  // TODO(nils): Figure out if we should we be aligning this as the
  // 20056 length SV instead, and if so, fix the aligner to do that.
  //
  // add_trace(64694328, 64694328 + 20054,
  //           "AAA");
  // run_vcf_test("4", "64694329", LongSequenceMatches("TATTGTGTATTATCACTAGTT.*"
  //                                                   "TATTATGTTATATTATATTATATAAGATATTATATTATATTTATGT"
  //                                                   "ATTAATATATTAATATGTTTATAT",
  //                                                   20054),
  //              "AAA", "1/1");

  // m_trace_enabled = true;
  // add_trace(64694330, 64694330 + 20052,
  //           "A");
  run_vcf_test("4", "64694331", LongSequenceMatches("TTGTGTATTATCACTAGTT.*"
                                                    "TATTATGTTATATTATATTATATAAGATATTATATTATATTTATGT"
                                                    "ATTAATATATTAATATGTTTATAT",
                                                    20052),
               "A", "1/1");

  return;
}

// This can cause a check failure if anchor lengths aren't normalized
// properly when doing anchor drop alignment.
TEST_F(hg002_test, sv_1_86356828) {
  m_options.min_overlap = 70;
  m_options.max_ploids = 50;

  call_at("1", "86356828", 150, 150);
}

// F0614 12:04:42.542244 49526 calc_coverage.cpp:130] Check failed:
// left_extent_end != right_extent_end (13219912 vs. 13219912)
// Assembly id=158755 [13052998:0+56112=13109110:0)
//
// This is an insert that happens between two extents:
//
//   "chr": "1",
//   "len": 9057730,
//   "name": "1:3995268",
//   "offset": 3995268

//   "chr": "1",
//   "len": 116914,
//   "name": "1:13102998",
//   "offset": 13102998

TEST_F(hg002_test, sv_1_13052998) {
  //
  call_region("1", 13052900, 13109200);
}

// 12 120875087 HG2_Ill_GATKHCSBGrefine_9070 T
// TTCAGGAGGCTGAGGCAGGGGAATCGCTTGAACCCGGGAGGCGGAGATTGCAGTGAGCTGAGATCGCGCCACTGCACTCCTCCAGCCTGGCAACAGAGCAAGATTCCGTTTCAAAAAAAAAAAAAAAAAAGTTTACTGTCATTTTATGTTATATACTTTTTTTTTAAAGTTTTATTTTTAAAGCTGCTTTTAGACAAGTCGAAGGAAGAAAGAAGGGGATAAGGAGGAAAGAATTTTGTAGACAAAATTTAACAGAGGTCAATTTTTTTTTTTTTTTTTTGTCTCCCAGGCTGGAGTGCAGTGGCATGATCTCAGCTCACTGCAACTTGCCCTCCCAGGTTCAAGCGATTCCTGTGCCTCAGCCACCCGGGTAGCTGGGATTACAGGTGTGCGCCACCATGCCCAGGTAATTTTTGTATTTTTAGTAGAGTTGGGGTTTCACCATGTTGGCCAGGCGGGTCTCGAACTTCTGACGTCAAGTGATCAGCCAGTCTCGGCATCCCAAAGTGCTGGGATTACAGGCGTGAACCACCACTCCCGGCCAGATGTCAATTTTTGTTTCCACAATTTCAAGGAAGAGAAAGCCAGTGTGACCAGAGGTCAAAAGATGAGAATGTTGGCCGGGTACGGTGGCTCATGCCTGTCATCCCACTATTTTGGGAGGCCGAGGCAGACAGATCACCTGAGGTCAAGAGTTTGAGACCAGCCTGGCTAACATGGTGAAACACCGTCTCTACTAAAAATACAAAAGAATTTGCTGGGCGTGGTGGTGCGTGCCTGTAATCCCAGCTAC
// 20 PASS
// END=120875087;SVTYPE=INS;SVLEN=792;ClusterIDs=HG2_PB_SVrefine2Falcon1plusDovetail_2680:HG4_Ill_SVrefine2DISCOVARDovetail_9705:HG2_Ill_SVrefine2DISCOVARplusDovetail_3034:HG4_Ill_GATKHCSBGrefine_9359:HG3_Ill_GATKHCSBGrefine_8960:HG2_Ill_GATKHCSBGrefine_9071:HG3_Ill_Spiral_14169:HG4_PB_HySA_19774:HG4_Ill_GATKHCSBGrefine_9358:HG3_Ill_GATKHCSBGrefine_8959:HG2_Ill_GATKHCSBGrefine_9070:HG2_PB_PB10Xdip_6594:HG2_PB_pbsv_13537:HG2_PB_HySA_24710:HG4_PB_pbsv_13536:HG3_PB_HySA_20450:HG3_PB_pbsv_13225:HG2_10X_SVrefine210Xhap12_9489:HG4_PB_SVrefine2PBcRDovetail_6759:HG3_Ill_SVrefine2DISCOVARDovetail_9838:HG2_PB_SVrefine2Falcon2Bionano_6724:HG3_PB_SVrefine2PBcRDovetail_6410:HG2_PB_SVrefine2PBcRplusDovetail_2549:HG2_PB_SVrefine2PB10Xhap12_9777:HG4_PB_SVrefine2Falcon1Dovetail_8250:HG3_PB_SVrefine2Falcon1Dovetail_8187;NumClusterSVs=26;ExactMatchIDs=HG2_Ill_GATKHCSBGrefine_9070:HG3_Ill_GATKHCSBGrefine_8959:HG4_Ill_GATKHCSBGrefine_9358;NumExactMatchSVs=3;ClusterMaxShiftDist=0.0665024630542;ClusterMaxSizeDiff=0.0665024630542;ClusterMaxEditDist=0.227465214761;PBcalls=15;Illcalls=10;TenXcalls=1;CGcalls=0;PBexactcalls=0;Illexactcalls=3;TenXexactcalls=0;CGexactcalls=0;HG2count=11;HG3count=8;HG4count=7;NumTechs=3;NumTechsExact=1;DistBack=-128;DistForward=4701;DistMin=-128;DistMinlt1000=TRUE;MultiTech=TRUE;MultiTechExact=FALSE;sizecat=300to999;DistPASSHG2gt49Minlt1000=FALSE;DistPASSMinlt1000=FALSE;MendelianError=FALSE;HG003_GT=1/1;HG004_GT=1/1;BREAKSIMLENGTH=255;REFWIDENED=12:120874912-120875166;REPTYPE=SIMPLEINS;TRall=FALSE;TRgt100=FALSE;TRgt10k=FALSE;segdup=FALSE;NumNeighbors=0;NumThresholdNeighbors=0
// GT:GTcons1:PB_GT:PB_REF:PB_ALT:PBHP_GT:PB_REF_HP1:PB_ALT_HP1:PB_REF_HP2:PB_ALT_HP2:TenX_GT:TenX_REF_HP1:TenX_ALT_HP1:TenX_REF_HP2:TenX_ALT_HP2:ILL250bp_GT:ILL250bp_REF:ILL250bp_ALT:ILLMP_GT:ILLMP_REF:ILLMP_ALT:BNG_LEN_DEL:BNG_LEN_INS:nabsys_svm
// 1/1:1/1:1/1:2:58:1/1:0:40:2:17:./.:1:8:0:3:./.:0:6:1/1:4:83:.:2315:.

TEST_F(hg002_test, HG2_Ill_GATKHCSBGrefine_9070) {
  call_region("12", 120875085, 120875088);
  EXPECT_THAT(
      m_assemblies,
      Contains(AllOf(
          VariantAt(120875086, 1,
                    "TTCAGGAGGCTGAGGCAGGGGAATCGCTTGAACCCGGGAGGCGGAGATTGCAGTGAGCTGAGATCGCGCCACTGCACT"
                    "CCTCCAGCCTGGCAACAGAGCAAGATTCCGTTTCAAAAAAAAAAAAAAAAAAGTTTACTGTCATTTTATGTTATATAC"
                    "TTTTTTTTTAAAGTTTTATTTTTAAAGCTGCTTTTAGACAAGTCGAAGGAAGAAAGAAGGGGATAAGGAGGAAAGAAT"
                    "TTTGTAGACAAAATTTAACAGAGGTCAATTTTTTTTTTTTTTTTTTGTCTCCCAGGCTGGAGTGCAGTGGCATGATCT"
                    "CAGCTCACTGCAACTTGCCCTCCCAGGTTCAAGCGATTCCTGTGCCTCAGCCACCCGGGTAGCTGGGATTACAGGTGT"
                    "GCGCCACCATGCCCAGGTAATTTTTGTATTTTTAGTAGAGTTGGGGTTTCACCATGTTGGCCAGGCGGGTCTCGAACT"
                    "TCTGACGTCAAGTGATCAGCCAGTCTCGGCATCCCAAAGTGCTGGGATTACAGGCGTGAACCACCACTCCCGGCCAGA"
                    "TGTCAATTTTTGTTTCCACAATTTCAAGGAAGAGAAAGCCAGTGTGACCAGAGGTCAAAAGATGAGAATGTTGGCCGG"
                    "GTACGGTGGCTCATGCCTGTCATCCCACTATTTTGGGAGGCCGAGGCAGACAGATCACCTGAGGTCAAGAGTTTGAGA"
                    "CCAGCCTGGCTAACATGGTGAAACACCGTCTCTACTAAAAATACAAAAGAATTTGCTGGGCGTGGTGGTGCGTGCCTG"
                    "TAATCCCAGCTAC" /* len=793 including anchor base */),
          /* GenotypeIs("1/1") */ _)))
      << dump_sv_assemblies(500);
}

// Per DEV-407, adam says this is a false positive that's causing problems:
// 10      4446063 .
// ATTACCCGGAAATGATCTGCTGAGACATAGGAAGCTATGCCCATCAGAGAAAACTAACAGAAGCCGGAATCTCCTCACATTCACGTTCACGACATCGCAGATATAGTCCGAGG
// A       100     PASS    NS=1;END=4446175;SVLEN=-112;SVTYPE=DEL
// GT:PG:GQ:PI:OV:DP:AD:PDP:PAD    0/1:0|1:43:1329:113:78:44,34:1:0,1
//
// TODO(nils): Make this test pass
TEST_F(hg002_test, DISABLED_fp_10_4446063) {
  //  add_assembly_trace(681);
  //  add_assembly_trace(1781);
  call_region("10", 4444400, 4447000);
  EXPECT_THAT(m_assemblies, Not(Contains(AllOf(VariantAt(4446062, 113, "A"), GenotypeIs("0/1")))));
  reset_assembly_trace();

  std::cout << "Ref before 4445996:\n";
  std::cout << get_ref_part_seq(4445996 - 100, 100) << "\n";
  std::cout << "Ref at 4445996->4446062, start of original traced assembly:\n";
  std::cout << get_ref_part_seq(4445996, 4446062 - 4445996) << "\n";
  std::cout << "Ref at 4446062->4446175, start of aligned assembly:\n";
  std::cout << get_ref_part_seq(4446062, 4446175 - 4446062) << "\n";
  std::cout << "Ref at 4446175->4446323, end of original traced assembly and aligned assembly:\n";
  std::cout << get_ref_part_seq(4446175, 4446323 - 4446175) << "\n";
  std::cout << "Ref after 4446623:\n";
  std::cout << get_ref_part_seq(4446323, 100) << "\n";
}

// 2 154752581       HG2_Ill_GATKHCSBGrefine_1548    A       AGGATGATTTTATATATATATATATATATATTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGGGATCTCGGCTCACTGCAAGCTCCGCCTCCCGGGTTAATGCCATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCGCCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCC       20      PASS    END=154752581;SVTYPE=INS;SVLEN=331;ClusterIDs=HG3_PB_HySA_3673:HG3_PB_HySA_3674:HG2_PB_pbsv_2376:HG3_PB_assemblyticsfalcon_14248:HG2_PB_assemblyticsPBcR_13820:HG2_Ill_breakscan11_3221:HG4_PB_pbsv_2387:HG2_PB_HySA_4512:HG3_PB_pbsv_2359:HG2_PB_PB10Xdip_16605:HG4_PB_HySA_3602:HG3_PB_SVrefine2Falcon1Dovetail_1466:HG3_PB_SVrefine2PBcRDovetail_1206:HG2_PB_SVrefine2Falcon2Bionano_964:HG2_PB_SVrefine2Falcon1plusDovetail_6405:HG2_PB_SVrefine2PB10Xhap12_1808:HG2_PB_SVrefine2PBcRplusDovetail_5932:HG2_10X_SVrefine210Xhap12_1807:HG4_Ill_MetaSV_2744:HG4_Ill_GATKHCSBGrefine_1610:HG2_Ill_MetaSV_2651:HG2_Ill_GATKHCSBGrefine_1548;NumClusterSVs=22;ExactMatchIDs=HG2_Ill_GATKHCSBGrefine_1548:HG2_Ill_MetaSV_2651:HG4_Ill_GATKHCSBGrefine_1610:HG4_Ill_MetaSV_2744;NumExactMatchSVs=4;ClusterMaxShiftDist=0.218890554723;ClusterMaxSizeDiff=0.218890554723;ClusterMaxEditDist=0.295031055901;PBcalls=16;Illcalls=5;TenXcalls=1;CGcalls=0;PBexactcalls=0;Illexactcalls=4;TenXexactcalls=0;CGexactcalls=0;HG2count=12;HG3count=6;HG4count=4;NumTechs=3;NumTechsExact=1;DistBack=35575;DistForward=-1;DistMin=-1;DistMinlt1000=TRUE;MultiTech=TRUE;MultiTechExact=FALSE;sizecat=300to999;DistPASSHG2gt49Minlt1000=FALSE;DistPASSMinlt1000=FALSE;MendelianError=FALSE;HG003_GT=1/1;HG004_GT=0/1;BREAKSIMLENGTH=14;REFWIDENED=2:154752582-154752595;REPTYPE=SIMPLEINS;TRall=FALSE;TRgt100=FALSE;TRgt10k=FALSE;segdup=FALSE;NumNeighbors=0;NumThresholdNeighbors=0       GT:GTcons1:PB_GT:PB_REF:PB_ALT:PBHP_GT:PB_REF_HP1:PB_ALT_HP1:PB_REF_HP2:PB_ALT_HP2:TenX_GT:TenX_REF_HP1:TenX_ALT_HP1:TenX_REF_HP2:TenX_ALT_HP2:ILL250bp_GT:ILL250bp_REF:ILL250bp_ALT:ILLMP_GT:ILLMP_REF:ILLMP_ALT:BNG_LEN_DEL:BNG_LEN_INS:nabsys_svm:PR:PA:PGT:PGQ 1/1:1/1:1/1:0:81:1/1:0:28:0:34:./.:0:5:0:2:1/1:0:17:1/1:3:75:.:.:.:3:11:1/1:16.09
// TODO(nils): Make this test pass]
TEST_F(hg002_test, DISABLED_manta_fn_HG2_Ill_GATKHCSBGrefine_1548) {
  // TODO(nils): Find out why this assembly isn't passing genotyping and fix it.
  // m_trace_enabled = true;
  m_options.min_read_depth = 0;
  m_options.min_depth_portion = 0;
  if (true) {
    // The two connections we care about making are:
    //
    // Read 1828669866, which starts with a lot of T's, and should end with read 1548396533 (not
    // 1070622507)

    // 1828669866: TTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGGGATCTCGGCTCAC
    // TGCAAGCTCCGCCTCCCGGGTTAATGCCATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCG ...
    // 1548396533: ... TGCAAGCTCCGCCTCCCGGGTTAATGCCATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCG
    // CCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATGGTCTCGAT

    // The other one is from 807920491 to 216970263 (not 1325438643):
    //
    // 807920491:
    // CGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCGGATGAT
    // TT ...
    // 216970263: ... ACACC
    // CGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCCGGAT

    // If we can make both of those connections in the same output assembly,
    // we should be able to output the proper variant.

    std::vector<uint32_t> trace_read_ids = {// Frist one:
                                            // First connection:
                                            1828669866, 1548396533,
                                            // Second connection:
                                            807920491, 216970263};

    if (true) {
      // Rev comp
      for (auto& read_id : trace_read_ids) {
        uint32_t rc_read_id = m_options.readmap->get_rev_comp(read_id);
        std::cout << "RC of read id " << read_id << " is: " << rc_read_id << "\n";
        read_id = rc_read_id;
      }
    }
    for (uint32_t read_id : trace_read_ids) {
      pop_tracer::add_debug_read(read_id);
    }
  }
  // g_trace_all_assemblies = true;
  // add_assembly_trace(199);
  // add_assembly_trace(200);
  // add_assembly_trace(201);
  // add_assembly_trace(202);
  // add_assembly_trace(203);

  // add_assembly_trace(670);
  // add_assembly_trace(671);
  // add_assembly_trace(672);
  // add_assembly_trace(673);
  // add_assembly_trace(674);

  // TODO(nils): Should genotype as "1/1", not "0/0":

  // Sequence including reference, for pileup, that we found:
  // TCTTGCAGTGTAATGATAAAAAAAGATATTGTCACACAAGAAATTCACAATTAAGTAAAAAGAGAATATAATCATAAAATAGTTTGTATGATATAATTATTTAAATAATAGTCTGCGAGTAAGAACAGGATAGTTGTTTCTGGAATTAAGGATGATTTTATATATATATATATATTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGGGATCTCGGCTCACTGCAAGCTCCGCCTCCCGGGTTAATGCCATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCGCCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCC
  // Sequence including refere,ce nfor pileup, that we should find:
  // TCTTGCAGTGTAATGATAAAAAAAGATATTGTCACACAAGAAATTCACAATTAAGTAAAAAGAGAATATAATCATAAAATAGTTTGTATGATATAATTATTTAAATAATAGTCTGCGAGTAAGAACAGGATAGTTGTTTCTGGAATTAAGGATGATTTTATATATATATATATATATATTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGGGATCTCGGCTCACTGCAAGCTCCGCCTCCCGGGTTAATGCCATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCGCCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCC

  add_print_seq_annotation("GGATGATTTTATA", "ref-repeat-in-variant");

  // It looks like:
  // Refereince is (ref before)-(ref-repeat-in-variant)-(ref-after)
  // Variant is (ref-before)-(ref-repeat-in-variant)-(more variant
  // stuff)-(ref-repeat-in-variant)-(ref-after)

  enable_annotated_sequences();

  add_print_seq_annotation("ATAGTTGTTTCTGGAATTAA", "ref-before");
  add_print_seq_annotation("TTTTTATGTTTTTTCCTAAA", "ref-after");
  add_print_seq_annotation(/* A here is wrong, G here is right */ "CCCGGCTAATTTTTTGTATTT",
                           "want-G-first");
  add_print_seq_annotation("TTTTATATATATATATATATATATTTTTT", "wrong-more-TA");

  run_vcf_test("2", "154752581", "A",
               "AGGATGATTTTATATATATATATATATATATTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGG"
               "AGTGCAGTGGCGGGATCTCGGCTCACTGCAAGCTCCGCCTCCCGGGTTAATGCCATTCTCCTGCCTCAGCCTCCCAAGTAGCT"
               "GGGACTACAGGCGCCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATG"
               "GTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGC"
               "C",
               "0/0");
  pop_tracer::clear_debug_reads();
}

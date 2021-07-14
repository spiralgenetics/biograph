#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>

#include "modules/variants/assemble_testutil.h"
#include "modules/variants/big_assemble_testutil.h"

using namespace testing;
using namespace variants;

constexpr bool k_full = false;

class hg001_downsample_test : public big_assemble_test {
 public:
  void SetUp() override {
    if (k_full) {
      // Non-downsampled version
      use_biograph("HG001_frc_8.bg");
    } else {
      // Downsampled version
      use_biograph("HG001.hs37d5.50x.11197.bg");
    }
  }
};

// This is a SV that wasn't found by bagpipe.
TEST_F(hg001_downsample_test, DISABLED_sv_20_50356899) {
  aoffset_t call_pos = 50357099;
  call_region("20", call_pos - 200, call_pos + 200 + 551);
  EXPECT_THAT(m_assemblies,
              Contains(AllOf(
                  VariantAt(call_pos - 1, 1,
                            "A"  // Vcf-style ref padding base
                            "ATTAGGCTGGGCACAGTGGCTCACACCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGTGGATCAC"
                            "CTGAGGTCGGGAGTTCAAGACCAGCCTAACCAACATGAGGAAACCCCGTCTCTACTAAAAATACAAAA"
                            "TTAGATGGGCGTGGTGGCGCATGCCTGTAATTCAAACTACTTGGAAGGCTGAGGCAGGAGAATTGCTT"
                            "GAACCCAGGAGACAGAGGTTGTGGTAAGCCAAGATCATGCCATTGTACTCCAGCATGGGCAACAAGAG"
                            "TGAGACTCCATCTCAAAAAAAAAAAAAATTAGCCAGGCGTGGTGGTGGGCACCTGTAATCCCAGCTAC"
                            "CCTGGAGACTGAGGCAGAAGAATCGCTTGAACCCAGGAGGCGGAGATTGCAGTGAGCCAAGATTACGC"
                            "CACTGCACTCCAGCCTGGGCACCAAGAGCAAAACCCTGTCTCAAAAAAATAAACAAATAAAAAGATTT"
                            "CTGTCTGCCACACGGCTGGGCCATGTGTAAAGACACATTCCTGTTGGTTTTATGTGTCTTGAATTCTA"
                            "ATGGG"),
                  GenotypeIs("1/1"))));
}

// Slow stats are:
// ambiguous_branch_cost: 188, ambiguous_ref_reads: 1, base_cost: 852,
// decrease_overlap_cost: 190, found_pairs: 488, increase_max_between_pair_cost:
// 36, matched_pairs: 33, max_assembly_len: 302, output_count: 5,
// pairs_used_cost: 33, search_not_fast_enough: 1, step_count: 210,
// too_many_ambiguous_steps: 1, too_many_pairs: 15
//
// TODO(nils): Figure out how to potentially speed up tracing around slow regions.
TEST_F(hg001_downsample_test, DISABLED_slow_region) {
  m_call_pos = 46400000;
  // call_region("16", 46400000, 46500000);
  call_region("16", 46400000, 46400010);
  EXPECT_LE(m_stats.step_count, 5) << m_stats;
}

// This can generate a check fail if the anchor sizes aren't handled right.
TEST_F(hg001_downsample_test, anchor_size_fail) { call_at("1", "4104619", 150, 150); }

// These are slow regions:
// (runtime 7 hours)
//  65.56 8 99600000 99700000
//  66.33 hs37d5 3600000 3700000
//  66.80 15 81900000 82000000
//  68.24 1 400000 500000
//  69.11 hs37d5 4900000 5000000
//  69.30 hs37d5 22800000 22900000
//  69.96 21 10500000 10600000
//  71.21 18 22200000 22300000
//  71.24 2 209600000 209700000
//  71.36 GL000214.1 0 100000
//  71.53 GL000224.1 0 100000
//  72.61 hs37d5 19300000 19400000
//  72.97 hs37d5 16600000 16700000
//  74.22 1 300000 400000
//  75.44 21 10400000 10500000
//  75.62 15 82000000 82100000
//  81.46 hs37d5 4400000 4500000
//  82.02 hs37d5 26500000 26600000
//  83.71 1 200000 300000
//  85.24 1 100000 200000
//  85.59 1 0 100000
//  85.88 hs37d5 2900000 3000000
//  88.48 NC_007605 0 100000
//  90.51 GL000218.1 0 100000
//  91.26 1 600000 700000
//  95.35 1 500000 600000
//  95.37 1 700000 800000
//  95.63 14 19000000 19100000
//  97.39 GL000199.1 0 100000
//  98.94 4 94500000 94600000
// 103.89 3 58800000 58900000
// 106.05 GL000199.1 100000 169874
// 106.28 X 61700000 61800000
// 111.28 MT 0 16569
// 121.69 GL000220.1 100000 161802
// 124.33 16 46300000 46400000
// 174.51 16 46400000 46500000
// 204.30 hs37d5 4700000 4800000
// 276.18 hs37d5 4800000 4900000
// 319.69 hs37d5 1400000 1500000

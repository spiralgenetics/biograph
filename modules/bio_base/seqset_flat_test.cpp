#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/io/file_io.h"
#include "modules/io/mem_io.h"
#include "modules/io/spiral_file_mem.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;
using namespace dna_testutil;

TEST(seqset_flat_test, seqset_flat_test) {
  std::vector<dna_sequence> test_seqs{tseq("abc"), tseq("bcd"), tseq("cde"),
                                      tseq("cdf"), tseq("dfg")};
  auto ss_f = seqset_for_reads(test_seqs);
  auto ss = &ss_f->get_seqset();

  // Copy down the set of all sequences from the original
  std::set<dna_sequence> expected_seqs;
  for (size_t i = 0; i < ss->size(); ++i) {
    std::cerr << "IDX " << i << ": "
              << ::testing::PrintToString(ss->ctx_entry(i).sequence()) << "\n";
    expected_seqs.insert(ss->ctx_entry(i).sequence());
  }

  spiral_file_mem_storage encoded;
  {
    seqset_flat_builder b(ss);
    spiral_file_create_mem c;
    b.build(c.create());
    encoded = c.close();
  }

  spiral_file_open_mem o(encoded);
  seqset_flat flat(o.open(), ss);

  std::set<dna_sequence> actual_seqs;
  for (dna_slice slice : flat) {
    actual_seqs.insert(dna_sequence(slice.begin(), slice.end()));
  }

  EXPECT_THAT(actual_seqs, ElementsAreArray(expected_seqs));
}

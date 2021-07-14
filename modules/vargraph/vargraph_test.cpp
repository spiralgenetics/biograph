
#include "vargraph.h"

#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/zip.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/readmap.h"

#include "modules/bio_base/seqset_testutil.h"
#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/reference_testutil.h"

#include <boost/algorithm/string.hpp>
#include <gtest/gtest.h>

using namespace dna_testutil;

typedef std::tuple<uint32_t, uint32_t, bool, dna_sequence, int, int, int, int> node_result;

std::string str_vec(const std::vector<int>& vec) {
  std::string r = "[";
  for(size_t i = 0; i < vec.size(); i++) {
    r += std::to_string(vec[i]);
    if (i != vec.size() - 1) { r += " "; }
  }
  r += "]";
  return r;
}

void dump_cov(const char* name, const vargraph::cov_info_t& ci) {
  SPLOG("  %s:", name);
  SPLOG("    base_cov: %s", str_vec(ci.base_cov).c_str());
  SPLOG("    span_cov: %s", str_vec(ci.span_cov).c_str());
};

void dump_graph(vargraph* vg, bool full_cov = true)
{
  SPLOG("Dump of vargraph with full_cov = %s", full_cov ? "true" : "false");
  for(const auto& kvp : vg->get_nodes()) {
    SPLOG("%s", kvp.second->as_string().c_str());
    if (full_cov) {
      dump_cov("unpaired", kvp.second->unpaired);
      dump_cov("paired", kvp.second->paired);
    }
  }
  for(const auto& e : vg->get_edges()) {
    SPLOG("%s->%s", e->upstream->as_string().c_str(), e->downstream->as_string().c_str());
    SPLOG("  unpaired: %u, paired: %u", e->unpaired, e->paired);
  }
}

void dump_reads(std::vector<std::vector<dna_sequence>> fake_reads) {
  for(auto const& read : fake_reads) {
    if (read.size() == 1) {
      SPLOG("unpr %s", read[0].as_string().c_str());
    } else if (read.size() == 2) {
      SPLOG("pair %s %s", read[0].as_string().c_str(), read[1].as_string().c_str());
    }
  }
}


vargraph::node_t* verify_node(vargraph* vg, uint32_t start, uint32_t end, bool is_ref, dna_sequence seq) {
  auto range = vg->get_nodes().equal_range(start);
  for(auto it = range.first; it != range.second; ++it) {
    vargraph::node_t* n = it->second;
    if (n->start == start &&
        n->end == end &&
        n->is_ref == is_ref &&
        n->seq == seq) {
      return n;
    }
  }
  EXPECT_TRUE(false) << "No matching node with sequence " << seq << " (is_ref=" << is_ref << ")";
  return nullptr;
}

void verify_edge_internal(vargraph* vg, vargraph::node_t* upstream, vargraph::node_t* downstream, uint32_t unpaired, uint32_t paired) {
  for(vargraph::edge_t* e : vg->get_edges()) {
    if (e->upstream == upstream && e->downstream == downstream) {
      EXPECT_EQ(e->unpaired, unpaired) << e->as_string();
      EXPECT_EQ(e->paired, paired) << e->as_string();
      return;
    }
  }
  EXPECT_TRUE(false) << "No matching node with upstream=" << upstream->as_string()
                     << " downstream=" << downstream->as_string();
}

#define VERIFY_EDGE(vg, upstream, downstream, unpaired, paired)                         \
  do {                                                                                  \
    SCOPED_TRACE("verify_edge(" #vg ", upstream=" #upstream ", downstream=" #downstream \
                 ", unpaired=" #unpaired ", paired=" #paired ")");                      \
    verify_edge_internal(vg, upstream, downstream, unpaired, paired);                   \
  } while (false)

TEST(vargraph, unpaired_hom_test) {
  /* Deletion >> abcd     jklm << Insertion
     @[50-60)     bcd-f    klmR   @[120, 120)
                   cd-fg    lmRn
                    d-fgh    mRno
                      fghi    Rnop
                               nopq
                 abcdefghijklmnopqrstuvwxyz << do a repl/snp also?
  */
  std::vector<std::vector<dna_sequence>> fake_reads =
                {
                  //unpaired
                  {tseq("abcd")},
                  {tseq("bcdf")},
                  {tseq("cdfg")},
                  {tseq("dfgh")},
                  {tseq("fghi")},
                  {tseq("jklm")},
                  {tseq("klmR")},
                  {tseq("lmRn")},
                  {tseq("mRno")},
                  {tseq("Rnop")},
                  {tseq("nopq")}
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijklmnopqrstuvwxyz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq);

  //Fake some variants
  //C insertion -- no coverage
  // When tsize > len - no error is thrown
  //vg->add_variant(24, 24, dna_C);
  //This doesn't work -- need to know what the dna_sequence expands to...
  //10bp deletion of 'c'
  vg->add_variant(40, 50, dna_sequence(""));
  //10bp insertion of R
  vg->add_variant(130, 130, tseq("R"));
  // Do trace
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  dump_graph(vg, true);

  // Verify nodes
  EXPECT_EQ(vg->get_nodes().size(), 5);
  auto n1 = verify_node(vg, 0, 40, true, tseq("abcd"));
  auto n2 = verify_node(vg, 40, 50, true, tseq("e"));
  auto n3 = verify_node(vg, 50, 130, true, tseq("fghijklm"));
  auto n4 = verify_node(vg, 130, 130, false, tseq("R"));
  auto n5 = verify_node(vg, 130, 260, true, tseq("nopqrstuvwxyz"));
  // Verify edges
  EXPECT_EQ(vg->get_edges().size(), 6);
  VERIFY_EDGE(vg, n1, n2, 0, 0);
  VERIFY_EDGE(vg, n2, n3, 0, 0);
  VERIFY_EDGE(vg, n1, n3, 3, 0); // Deletion ede
  VERIFY_EDGE(vg, n3, n4, 3, 0);
  VERIFY_EDGE(vg, n4, n5, 3, 0);
  VERIFY_EDGE(vg, n3, n5, 0, 0);
  // Verify some interesting coverage points
  // Check that span coverage drops to 0 at the exact point where fghi and ijkl fail to overlap
  EXPECT_EQ(n3->unpaired.span_cov[38], 1);
  EXPECT_EQ(n3->unpaired.span_cov[39], 0);
  EXPECT_EQ(n3->unpaired.span_cov[40], 1);
  // Check base coverage on inserted base (which is only good if reads split across nodes)
  EXPECT_EQ(n4->unpaired.base_cov[5], 4);
  // Check deleted node coverage is 0
  EXPECT_EQ(n2->unpaired.base_cov[5], 0);

  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}


// TEST(vargraph, paired_hom_test) {
//   /* Deletion >> abcd     jklm << Insertion
//      @[50-60)     bcd-f    klmR   @[120, 120)
//                    cd-fg    lmRn
//                     d-fgh    mRno
//                       fghi    Rnop
//                                nopq    wxyz < pair for all
//                                               so a:i don't have
//                  abcdefghijklmnopqrstuvwxyz
//                           | mx 200 insert |
//   */

//   std::vector<std::vector<dna_sequence>> fake_reads =
//                 {
//                   //unpaired -- Need to check directionality, too
//                   //paired - outof distance
//                   {tseq("abcd"), tseq_rc("wxyz")},
//                   {tseq("bcdf"), tseq_rc("wxyz")},
//                   {tseq("cdfg"), tseq_rc("wxyz")},
//                   {tseq("dfgh"), tseq_rc("wxyz")},
//                   {tseq("fghi"), tseq_rc("wxyz")},
//                   //paired - in distance
//                   {tseq("jklm"), tseq_rc("wxyz")},
//                   {tseq("klmR"), tseq_rc("wxyz")},
//                   {tseq("lmRn"), tseq_rc("wxyz")},
//                   {tseq("mRno"), tseq_rc("wxyz")},
//                   {tseq("Rnop"), tseq_rc("wxyz")},
//                   {tseq("nopq"), tseq_rc("wxyz")},
//                  };

//   //Fake reference & BioGraph
//   std::unique_ptr<const reference> ref;
//   dna_sequence ref_seq = tseq("abcdefghijklmnopqrstuvwxyz");
//   int tsize = ref_seq.size();
//   ref = create_reference({ref_seq});

//   auto m_biograph = biograph_for_reads(fake_reads);
//   std::shared_ptr<const seqset> m_seqset;
//   m_seqset = std::move(m_biograph.first);
//   auto m_readmap = std::move(m_biograph.second);

//   // Build VarGraph
//   vargraph *vg = nullptr;
//   vg = new vargraph(ref_seq, 10, 200);

//   //Fake some variants
//   //C insertion -- no coverage
//   // When tsize > len - no error is thrown
//   //vg->add_variant(24, 24, dna_C);
//   //This doesn't work -- need to know what the dna_sequence expands to...
//   //10bp deletion of 'c'
//   vg->add_variant(40, 50, dna_sequence("")); // deletion doesn't make a 'node' anymore
//   //10bp insertion of R
//   vg->add_variant(130, 130, tseq("R"));
//   //Then trace
//   vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

//   //Verify nodes exist
//   EXPECT_EQ(vg->get_nodes().size(), 5);
//   auto n1 = verify_node(vg, 0, 40, true, tseq("abcd"));
//   auto n2 = verify_node(vg, 40, 50, true, tseq("e"));
//   auto n3 = verify_node(vg, 50, 130, true, tseq("fghijklm"));
//   auto n4 = verify_node(vg, 130, 130, false, tseq("R"));
//   auto n5 = verify_node(vg, 130, 260, true, tseq("nopqrstuvwxyz"));

//   // Verify edges
//   EXPECT_EQ(vg->get_edges().size(), 6);
//   VERIFY_EDGE(vg, n1, n2, 0, 0);
//   VERIFY_EDGE(vg, n2, n3, 0, 0);
//   VERIFY_EDGE(vg, n1, n3, 3, 0); // Deletion edge - all unpaired
//   VERIFY_EDGE(vg, n3, n4, 0, 3); // Insertion edge in - all paired
//   VERIFY_EDGE(vg, n4, n5, 0, 3); // Insertion edge out - all paired
//   VERIFY_EDGE(vg, n3, n5, 0, 0); // Ref - no support

//   //should put interesting point checks
//   dump_graph(vg);
//   EXPECT_EQ(n3->unpaired.span_cov[0], 4);
//   EXPECT_EQ(n3->paired.span_cov[0], 0);
//   EXPECT_EQ(n4->paired.span_cov[0], 4);
//   EXPECT_EQ(n4->paired.span_cov[8], 4);
//   EXPECT_EQ(n4->unpaired.span_cov[9], 4);
//   EXPECT_EQ(n5->paired.base_cov[0], 4);
//   EXPECT_EQ(n5->unpaired.base_cov[0], 0);
//   EXPECT_EQ(n5->unpaired.span_cov[128], 5);
//   EXPECT_EQ(n5->paired.span_cov[128], 6);

//   //Then look at the results
//   dump_graph(vg);
//   dump_reads(fake_reads);
//   SPLOG("ref %s", ref_seq.as_string().c_str());

//   delete vg;
// }

TEST(vargraph, redundant_read_test) {
  /* what happens if we have redundancy in the ref sequence?
     on one copy, we'll have proper pairs, on another we'll have improper pairs
     Does the vargraph 1) only use the read once? if so where?
        2) use the read is used twice, once proper, once improper, are they correctly paired/unpaired

  ref abcdefghijklcdefgmnopqrstuvwxyz
        ^rep1     ^rep2

                 lcde   nopq
                  cdef   opqr
                   defg   pqrs
                    efgm   qrst
        Max Insert 150
        Min Insert 50 should be fine
  */
  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijklcdefgmnopqrstuvwxy");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  std::vector<std::vector<dna_sequence>> fake_reads =
                {
                  {tseq("lcde"), tseq_rc("nopq")},
                  {tseq("cdef"), tseq_rc("opqr")}, //r1 repeat
                  {tseq("defg"), tseq_rc("pqrs")}, //r1 repeat
                  {tseq("efgm"), tseq_rc("qrst")}
                };

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 50, 150);
  // No variants
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  dump_graph(vg);
  EXPECT_EQ(vg->get_nodes().size(), 1);
  auto n1 = verify_node(vg, 0, 300, true, ref_seq);
  //The repeats make the cdef used twice
  EXPECT_EQ(n1->paired.base_cov[30], 0);
  EXPECT_EQ(n1->unpaired.base_cov[30], 2); // the l read's d not counted
  EXPECT_EQ(n1->paired.base_cov[130], 3);
  EXPECT_EQ(n1->unpaired.base_cov[130], 0);
}

TEST(vargraph, het_test) {
  /* Deletion >> abcd     jklm << Insertion
     @[50-60)     bcd-f    klmR   @[120, 120)
                   cd-fg    lmRn
                    d-fgh    mRno
                      fghi    Rnop
                               nopq    wxyz < pair for all up to now
                                              so a:i don't have proper pair
                  bcde      klmn
                   cdef      lmno
                    defg      mnop  < pairs supporting refernce
                     efgh    lmno   <  unpaired supporting reference
                 abcdefghijklmnopqrstuvwxyz
                          | mx 200 insert |

  */

  std::vector<std::vector<dna_sequence>> fake_reads =
                {
                  //paired - out of distance
                  {tseq("abcd"), tseq_rc("wxyz")},
                  {tseq("bcdf"), tseq_rc("wxyz")},
                  {tseq("cdfg"), tseq_rc("wxyz")},
                  {tseq("dfgh"), tseq_rc("wxyz")},
                  {tseq("fghi"), tseq_rc("wxyz")},
                  //paired - in distance
                  {tseq("jklm"), tseq_rc("wxyz")},
                  {tseq("klmR"), tseq_rc("wxyz")},
                  {tseq("lmRn"), tseq_rc("wxyz")},
                  {tseq("mRno"), tseq_rc("wxyz")},
                  {tseq("Rnop"), tseq_rc("wxyz")},
                  {tseq("nopq"), tseq_rc("wxyz")},
                  //paired ref allele reads
                  {tseq("bcde"), tseq_rc("klmn")},
                  {tseq("cdef"), tseq_rc("lmno")},
                  {tseq("defg"), tseq_rc("mnop")},
                  //unpaired ref allele reads
                  {tseq("efgh")},
                  {tseq_rc("lmno")},
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijklmnopqrstuvwxyz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 10, 200);

  //Fake some variants
  //C insertion -- no coverage
  // When tsize > len - no error is thrown
  //vg->add_variant(24, 24, dna_C);
  //This doesn't work -- need to know what the dna_sequence expands to...
  //10bp deletion of 'c'
  vg->add_variant(40, 50, dna_sequence("")); // deletion doesn't make a 'node' anymore
  //10bp insertion of R
  vg->add_variant(130, 130, tseq("R"));
  //Then trace
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  //Verify nodes exist
  ASSERT_EQ(vg->get_nodes().size(), 5);
  auto n1 = verify_node(vg, 0, 40, true, tseq("abcd"));
  auto n2 = verify_node(vg, 40, 50, true, tseq("e"));
  auto n3 = verify_node(vg, 50, 130, true, tseq("fghijklm"));
  auto n4 = verify_node(vg, 130, 130, false, tseq("R"));
  auto n5 = verify_node(vg, 130, 260, true, tseq("nopqrstuvwxyz"));

  // Verify edges
  ASSERT_EQ(vg->get_edges().size(), 6);
  VERIFY_EDGE(vg, n1, n2, 0, 3);
  VERIFY_EDGE(vg, n2, n3, 1, 2); // Ref allele
  VERIFY_EDGE(vg, n1, n3, 3, 0); // Deletion edge - all unpaired
  VERIFY_EDGE(vg, n3, n4, 0, 3); // Insertion edge in - all paired
  VERIFY_EDGE(vg, n4, n5, 0, 3); // Insertion edge out - all paired
  VERIFY_EDGE(vg, n3, n5, 1, 3); // Ref allele

  //should put interesting point checks
  dump_graph(vg);
  EXPECT_EQ(n3->unpaired.span_cov[0], 5);
  EXPECT_EQ(n3->paired.span_cov[0], 2);
  EXPECT_EQ(n4->paired.span_cov[0], 4);
  EXPECT_EQ(n4->paired.span_cov[8], 4);
  EXPECT_EQ(n5->paired.base_cov[0], 7);
  EXPECT_EQ(n5->unpaired.base_cov[0], 1);
  EXPECT_EQ(n5->unpaired.span_cov[128], 5);
  EXPECT_EQ(n5->paired.span_cov[128], 6);

  //Then look at the results
  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}

TEST(vargraph, multi_pair) {
  /*
  What happens if a read maps twice within a node?


  abcdefghijklmnopqrstjklmnz
   bcd-f   jklm       jklm
   readaln1  read2aln    read2aln2
  */

  std::vector<std::vector<dna_sequence>> fake_reads =
                {
                  //paired - both distance
                  {tseq("bcdf"), tseq_rc("jklm")},
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijklmnopqrstjklmnz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 10, 200);

  //Fake some variants
  //C insertion -- no coverage
  // When tsize > len - no error is thrown
  //vg->add_variant(24, 24, dna_C);
  //This doesn't work -- need to know what the dna_sequence expands to...
  //10bp deletion of 'c'
  vg->add_variant(40, 50, dna_sequence("")); // deletion doesn't make a 'node' anymore

  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  //Verify nodes exist
  ASSERT_EQ(vg->get_nodes().size(), 3);
  auto n1 = verify_node(vg, 0, 40, true, tseq("abcd"));
  auto n2 = verify_node(vg, 40, 50, true, tseq("e"));
  auto n3 = verify_node(vg, 50, 260, true, tseq("fghijklmnopqrstjklmnz"));

  ASSERT_EQ(vg->get_edges().size(), 3);
  VERIFY_EDGE(vg, n1, n2, 0, 0);
  VERIFY_EDGE(vg, n2, n3, 0, 0);
  VERIFY_EDGE(vg, n1, n3, 0, 1);

  //should put interesting point checks
  dump_graph(vg);

  //Then look at the results
  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}

TEST(vargraph, insert_pair_test) {
  /*
  Can I break the insert?
             /G\
  abcdefghijkl  mnopqrstuvwxyz
             \A/
           jklGm       tuvw   : fail - insertsize is 151
           jklGmnopq          : fail - 81
           jklAm   pqrs       : pass 101

   min_pair_dist = 101
   max_pair_dist = 101
   bcde should fail. pqrs should pass. nopq should fail

   readaln1  read2aln    read2aln2
  */

  std::vector<std::vector<dna_sequence>> fake_reads =
                { //paired - both distance
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("tuvw"),  },
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("nopq"),  },
                  {tseq("jkl") + dna_A + tseq("m"), tseq_rc("pqrs") },
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijkl") + dna_G + tseq("mnopqrstuvwxyz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 101, 101);

  //Fake some variants - not needed since it's not indel. het snp
  vg->add_variant(120, 121, dna_A);
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  //Verify nodes exist
  ASSERT_EQ(vg->get_nodes().size(), 4);
  auto n1 = verify_node(vg, 0, 120, true, tseq("abcdefghijkl"));
  auto n2 = verify_node(vg, 120, 121, true, dna_G);
  auto n3 = verify_node(vg, 120, 121, false, dna_A);
  auto n4 = verify_node(vg, 121, 261, true, tseq("mnopqrstuvwxyz"));

  ASSERT_EQ(vg->get_edges().size(), 4);
  VERIFY_EDGE(vg, n1, n2, 2, 0);
  VERIFY_EDGE(vg, n1, n3, 0, 1);
  VERIFY_EDGE(vg, n2, n4, 2, 0);
  VERIFY_EDGE(vg, n3, n4, 0, 1);

  //Then look at the results
  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}


TEST(vargraph, insert_pair_test2) {
  /*
  Can I break the insert?
              /G\  /-----\
  abcdefghijkl   mn-opqrs-tuvwxyz
              \A/
 read1     jkl G m         uvwx  : pass - insertsize is 101 because I'll make a 50bp del
 read2     jkl G mn opq          : fail - 81
 read3     jkl A m   pqrs        : pass 101

   min_pair_dist = 101
   max_pair_dist = 101

  OKAY! When I run this with 151 max (insert of read1 going through reference allele over deletion)
        I get correct paired coverage on the G ref allele
  BUT I do not get correct paired coverage on it's up/dn node - so the coverage goes 2-1-2!!

  ALSO! When I run this test with 101 max (insert of read1 going through the deletion allele) it stays 2-2-2.

  SO! this is where the problem is. It's not properly marking pairs on all reads in all nodes!
  I'm suspicious of the reachable
  */

  std::vector<std::vector<dna_sequence>> fake_reads =
                { //paired - both distance
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("uvwx"),  },
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("nopq"),  },
                  {tseq("jkl") + dna_A + tseq("m"), tseq_rc("pqrs") },
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijkl") + dna_G + tseq("mnopqrstuvwxyz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 101, 200);

  //Fake some variants - not needed since it's not indel. het snp
  vg->add_variant(120, 121, dna_A);
  vg->add_variant(141, 191, dna_sequence(""));
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  //Verify nodes exist
  ASSERT_EQ(vg->get_nodes().size(), 6);
  auto jkl = verify_node(vg, 0, 120, true, tseq("abcdefghijkl"));
  auto g_insert = verify_node(vg, 120, 121, true, dna_G);
  auto a_insert = verify_node(vg, 120, 121, false, dna_A);
  auto mn = verify_node(vg, 121, 141, true, tseq("mn"));
  auto opqrs = verify_node(vg, 141, 191, true, tseq("opqrs"));
  auto tuv = verify_node(vg, 191, 261, true, tseq("tuvwxyz"));

  ASSERT_EQ(vg->get_edges().size(), 7);
  VERIFY_EDGE(vg, jkl, g_insert, 1, 1);
  VERIFY_EDGE(vg, jkl, a_insert, 0, 1);
  VERIFY_EDGE(vg, g_insert, mn, 1, 1);
  VERIFY_EDGE(vg, a_insert, mn, 0, 1);
  VERIFY_EDGE(vg, mn, opqrs, 1, 0);
  VERIFY_EDGE(vg, opqrs, tuv, 0, 0);

  //Then look at the results
  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}

TEST(vargraph, insert_pair_test3) {
  /*
  This is the same as insert_pair_test2 except I'm going to have read1's dn pair go over the del edge
  It also fails for the same reason. So hitting a new edge doesn't update the info
              /G\  /-----\
  abcdefghijkl   mn-opqrs-tuvwxyz
              \A/
 read1     jkl G m      s tuv    : pass - insertsize is 101 because I'll make a 50bp del
 read2     jkl G mn opq          : fail - 81
 read3     jkl A m   opqrs       : pass 101

   min_pair_dist = 101
   max_pair_dist = 101

  OKAY! When I run this with 151 max (insert of read1 going through reference allele over deletion)
        I get correct paired coverage on the G ref allele
  BUT I do not get correct paired coverage on it's up/dn node - so the coverage goes 2-1-2!!

  ALSO! When I run this test with 101 max (insert of read1 going through the deletion allele) it stays 2-2-2.

  SO! this is where the problem is. It's not properly marking pairs on all reads in all nodes!
  I'm suspicious of the reachable
  */

  std::vector<std::vector<dna_sequence>> fake_reads =
                { //paired - both distance
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("stuv"),  },
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("nopq"),  },
                  {tseq("jkl") + dna_A + tseq("m"), tseq_rc("pqrs") },
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijkl") + dna_G + tseq("mnopqrstuvwxyz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 101, 200);

  //Fake some variants - not needed since it's not indel. het snp
  vg->add_variant(120, 121, dna_A);
  vg->add_variant(141, 191, dna_sequence(""));
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  //Verify nodes exist
  ASSERT_EQ(vg->get_nodes().size(), 6);
  auto n1 = verify_node(vg, 0, 120, true, tseq("abcdefghijkl"));
  auto n2 = verify_node(vg, 120, 121, true, dna_G);
  auto n3 = verify_node(vg, 120, 121, false, dna_A);
  auto n4 = verify_node(vg, 121, 141, true, tseq("mn"));
  auto n5 = verify_node(vg, 141, 191, true, tseq("opqrs"));
  auto n6 = verify_node(vg, 191, 261, true, tseq("tuvwxyz"));

  ASSERT_EQ(vg->get_edges().size(), 7);
  VERIFY_EDGE(vg, n1, n2, 1, 1);
  VERIFY_EDGE(vg, n1, n3, 0, 1);
  VERIFY_EDGE(vg, n2, n4, 1, 1);
  VERIFY_EDGE(vg, n3, n4, 0, 1);
  VERIFY_EDGE(vg, n4, n5, 1, 0);
  VERIFY_EDGE(vg, n5, n6, 0, 1);

  //Then look at the results
  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}

TEST(vargraph, insert_culprit) {
  /*
  Just the bad read - I'm going to dump everything from vg to see what it does
  Dude -- this one passes... so now I need to check if I can make the above pass with
  non-identical reads..
              /G\  /-----\
  abcdefghijkl   mn-opqrs-tuvwxyz
              \A/
 read1     jkl G m         uvwx  : pass - insertsize is 151 because I'll make a 50bp del

   min_pair_dist = 101
   max_pair_dist = 101
  */

  std::vector<std::vector<dna_sequence>> fake_reads =
                { //paired - both distance
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("uvwx") },
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijkl") + dna_G + tseq("mnopqrstuvwxyz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 101, 200);

  //Fake some variants - not needed since it's not indel. het snp
  vg->add_variant(120, 121, dna_A);
  vg->add_variant(141, 191, dna_sequence(""));
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  //Verify nodes exist
  ASSERT_EQ(vg->get_nodes().size(), 6);
  auto n1 = verify_node(vg, 0, 120, true, tseq("abcdefghijkl"));
  auto n2 = verify_node(vg, 120, 121, true, dna_G);
  auto n3 = verify_node(vg, 120, 121, false, dna_A);
  auto n4 = verify_node(vg, 121, 141, true, tseq("mn"));
  auto n5 = verify_node(vg, 141, 191, true, tseq("opqrs"));
  auto n6 = verify_node(vg, 191, 261, true, tseq("tuvwxyz"));

  ASSERT_EQ(vg->get_edges().size(), 7);
  VERIFY_EDGE(vg, n1, n2, 0, 1);
  VERIFY_EDGE(vg, n1, n3, 0, 0);
  VERIFY_EDGE(vg, n2, n4, 0, 1);
  VERIFY_EDGE(vg, n3, n4, 0, 0);
  VERIFY_EDGE(vg, n4, n5, 0, 0);
  VERIFY_EDGE(vg, n5, n6, 0, 0);

  //Then look at the results
  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}


TEST(vargraph, insert_culprit_alibi) {
  /*
  This is the same as insert_culprit and insert_test_2, but I'm
  changing the jkl read so that it's mapping to a unique read.
  If this passes, I'm effed
              /G\  /-----\
  abcdefghijkl   mn-opqrs-tuvwxyz
              \A/
 read1     jkl G m         uvwx  : pass - insertsize is 151 or min 101
 read3      kl A mn   qrst       : pass - 101

   min_pair_dist = 101
   max_pair_dist = 200

  OKAY! When I run this with 151 max (insert of read1 going through reference allele over deletion)
        I get correct paired coverage on the G ref allele
  BUT I do not get correct paired coverage on it's up/dn node - so the coverage goes 2-1-2!!

  ALSO! When I run this test with 101 max (insert of read1 going through the deletion allele) it stays 2-2-2.

  SO! this is where the problem is. It's not properly marking pairs on all reads in all nodes!
  I'm suspicious of the reachable
  */

  std::vector<std::vector<dna_sequence>> fake_reads =
                { //paired - both distance
                  {tseq("jkl") + dna_G + tseq("m"), tseq_rc("uvwx"),  },
                  {tseq("kl") + dna_A + tseq("mn"), tseq_rc("qrst") },
                 };

  //Fake reference & BioGraph
  std::unique_ptr<const reference> ref;
  dna_sequence ref_seq = tseq("abcdefghijkl") + dna_G + tseq("mnopqrstuvwxyz");
  int tsize = ref_seq.size();
  ref = create_reference({ref_seq});

  auto m_biograph = biograph_for_reads(fake_reads);
  std::shared_ptr<const seqset> m_seqset;
  m_seqset = std::move(m_biograph.first);
  auto m_readmap = std::move(m_biograph.second);

  // Build VarGraph
  vargraph *vg = nullptr;
  vg = new vargraph(ref_seq, 101, 101);

  //Fake some variants - not needed since it's not indel. het snp
  vg->add_variant(120, 121, dna_A);
  vg->add_variant(141, 191, dna_sequence(""));
  vg->trace(*m_seqset, *m_readmap.get(), 0, tsize);

  //Verify nodes exist
  ASSERT_EQ(vg->get_nodes().size(), 6);
  auto n1 = verify_node(vg, 0, 120, true, tseq("abcdefghijkl"));
  auto n2 = verify_node(vg, 120, 121, true, dna_G);
  auto n3 = verify_node(vg, 120, 121, false, dna_A);
  auto n4 = verify_node(vg, 121, 141, true, tseq("mn"));
  auto n5 = verify_node(vg, 141, 191, true, tseq("opqrs"));
  auto n6 = verify_node(vg, 191, 261, true, tseq("tuvwxyz"));

  ASSERT_EQ(vg->get_edges().size(), 7);
  VERIFY_EDGE(vg, n1, n2, 0, 1);
  VERIFY_EDGE(vg, n1, n3, 0, 1);
  VERIFY_EDGE(vg, n2, n4, 0, 1);
  VERIFY_EDGE(vg, n3, n4, 0, 1);
  VERIFY_EDGE(vg, n4, n5, 0, 0);
  VERIFY_EDGE(vg, n5, n6, 0, 1);

  //Then look at the results
  dump_graph(vg);
  dump_reads(fake_reads);
  SPLOG("ref %s", ref_seq.as_string().c_str());

  delete vg;
}

#include "modules/test/test_utils.h"
#include "modules/bio_base/call_structural.h"
#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/bio_base/reference.h"
#include <gtest/gtest.h>

TEST(left_align, a)
{
	compound_cost cost1(10, 0);
	compound_cost cost2(10, 20);
	compound_cost cost3(11, 0);

	ASSERT_TRUE( cost1 < cost2 );
	ASSERT_TRUE( cost2 < cost3 );
	ASSERT_TRUE( cost2 == compound_cost(10,20) );

	compound_cost cost4;
	compound_cost cost5 = cost4 + cost2;

	ASSERT_TRUE( cost5 == cost2 );
}

TEST(call_structural, left_align_2b_del) 
{
	//                012345678
	dna_sequence ref("ACTCACATCCA");
	dna_sequence seq("ACTCATCCA");
	//                012345678

	std::vector<sv_out> out = call_structural(seq, dna_slice(ref), dna_slice(ref), compound_cost(100));

	ASSERT_FALSE( out[0].is_structural );

	ASSERT_EQ( out[0].seq_begin,	3 );				// before left-alignment: 5
	ASSERT_EQ( out[0].seq_end,		3 );				// before left-alignment: 5
	ASSERT_EQ( out[0].left_ref ,	ref.begin() + 2 );	// before left-alignment: 4
	ASSERT_EQ( out[0].right_ref ,	ref.begin() + 5 );	// before left-alignment: 7
}

TEST(call_structural, left_align_reverse_complement) 
{
	//                012345678
	dna_sequence ref("ACTCACATCCA");
	dna_sequence refc = ref.rev_comp();
	// we are deleting 2 base pairs CACA -> CA and aligning against the reverse complement (complement() returns the reverse comp!!)
	// So CACA gets turns into TGTG
	// Since we want to left-align the deletion, that means the first TG at position 6-5 should go.
	//                TGGATGTGAGT
	//                 9876543210
	dna_sequence seq("ACTCATCCA");
	//                012345678

	dna_slice refc_c = dna_slice(refc).rev_comp();
	std::vector<sv_out> out = call_structural(seq, refc_c, refc_c, compound_cost(100));

	ASSERT_FALSE( out[0].is_structural );

	ASSERT_EQ( out[0].seq_begin,	5 );
	ASSERT_EQ( out[0].seq_end,		5 );
	//SPLOG("left_ref - refc_c.begin = %ld", out[0].left_ref - refc_c.begin);
	//SPLOG("right_ref - refc_c.begin = %ld", out[0].right_ref - refc_c.begin);
	ASSERT_EQ( out[0].left_ref ,	refc_c.begin() + 4 );
	ASSERT_EQ( out[0].right_ref ,	refc_c.begin() + 7 );
}

TEST(call_structural, left_align_4b_ins) 
{
	//                012345678
	dna_sequence ref("ACTCATCCA");
	dna_sequence seq("ACTCACACATCCA");
	//                0123456789

	std::vector<sv_out> out = call_structural(seq, dna_slice(ref), dna_slice(ref), compound_cost(100));

	ASSERT_FALSE( out[0].is_structural );

	ASSERT_EQ( out[0].seq_begin,	3 );
	ASSERT_EQ( out[0].seq_end,		7 );
	ASSERT_EQ( out[0].left_ref ,	ref.begin() + 2 );
	ASSERT_EQ( out[0].right_ref ,	ref.begin() + 3 );
}

TEST(call_structural, test1) 
{
	//                012345678
	dna_sequence ref("ACTCATCCA");
	dna_sequence seq("ACTTTCCA");
	//                01234567
	std::vector<sv_out> out = call_structural(seq, dna_slice(ref), dna_slice(ref), compound_cost(100));
	ASSERT_EQ(out.size(), size_t(1));
	ASSERT_TRUE(out[0].is_structural == false);
	ASSERT_TRUE(out[0].seq_begin == 3);
	ASSERT_TRUE(out[0].seq_end == 4);
	ASSERT_TRUE(out[0].left_ref == ref.begin() + 2);
	ASSERT_TRUE(out[0].right_ref == ref.begin() + 5);
}

TEST(call_structural, test2) 
{
	dna_sequence ref1("GTGACTCATCCAAC");
        dna_sequence ref2("TCCCACAAGCACGG");
	dna_sequence seq ("GTGACTTGGCACGG");
	//                 0123456789
	std::vector<sv_out> out = call_structural(seq, dna_slice(ref1), dna_slice(ref2), compound_cost(100));
	ASSERT_EQ(out.size(), size_t(1));
	ASSERT_TRUE(out[0].is_structural == true);
	ASSERT_TRUE(out[0].seq_begin == 6);
	ASSERT_TRUE(out[0].seq_end == 8);
	ASSERT_TRUE(out[0].left_ref == ref1.begin() + 5);
	ASSERT_TRUE(out[0].right_ref == ref2.begin() + 8);
}

TEST(call_structural, test3) 
{
	//                 01234567890123456789012345678901234567890123
	dna_sequence ref1("GACTCACGCTGACTAGCTACCCCATACCAGGGGCATATCATCCA");
	//                 ====*===========***
	dna_sequence seq ("GACTTACGCTGACTAGGGGGCCCATTACGCGGCAGTCCAACAC");
        //                                 ***===============**.\\\\\\\  //
        dna_sequence ref2("CCACAAGCCACAGGATCGAGCCCATTACGCGGCATACCCAACAC");
	//                 01234567890123456789012345678901234567890123
	std::vector<sv_out> out = call_structural(seq, dna_slice(ref1), dna_slice(ref2), compound_cost(100));
	ASSERT_EQ(out.size(), size_t(3));
	ASSERT_TRUE(out[0].is_structural == false);
	ASSERT_TRUE(out[0].seq_begin == 4);
	ASSERT_TRUE(out[0].seq_end == 5);
	ASSERT_TRUE(out[0].left_ref == ref1.begin() + 3);
	ASSERT_TRUE(out[0].right_ref == ref1.begin() + 5);
	ASSERT_TRUE(out[1].is_structural == true);
	ASSERT_TRUE(out[1].seq_begin == 16);
	ASSERT_TRUE(out[1].seq_end == 19);
	ASSERT_TRUE(out[1].left_ref == ref1.begin() + 15);
	ASSERT_TRUE(out[1].right_ref == ref2.begin() + 19);
	ASSERT_TRUE(out[2].is_structural == false);
	ASSERT_TRUE(out[2].seq_begin == 34);
	ASSERT_TRUE(out[2].seq_end == 36);
	ASSERT_TRUE(out[2].left_ref == ref2.begin() + 33);
	ASSERT_TRUE(out[2].right_ref == ref2.begin() + 37);
}

TEST(call_structural, test4)
{
	
	dna_sequence ref1("GACTCACGCTGACTAGCTACCAGGGGCATATCATCCA");
	//                 0123456789012345***********
	dna_sequence seq( "GACTCACGCTGACTAGTTCAGGCATCATCTTAAGTCTGACTCGATTTTACACCCTCGTTGAATACCCCATACCAGGGGCATATCATCCA");
	//                           1         2         3         4         5         6         7
	//                 012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
	//                                                                                 *****|
	dna_sequence ref2("GACTCACGCTGACTAGCTACCAGGGGCATATCATCCA");
	//                 012345678901234567
	std::vector<sv_out> out = call_structural(seq, dna_slice(ref1), dna_slice(ref2), 10);
	ASSERT_EQ(out.size(), size_t(1));
	ASSERT_TRUE(out[0].is_structural == true);
	ASSERT_TRUE(out[0].seq_begin == 16);
	ASSERT_TRUE(out[0].seq_end == 69);
	ASSERT_TRUE(out[0].left_ref == ref1.begin() + 15);
	ASSERT_TRUE(out[0].right_ref == ref2.begin() + 17);
}

TEST(call_structural, DISABLED_evil)
{
	reference ref("homo_sapiens_GRCh37");
	dna_sequence assembled("TTTTGCTATGCAGAAGCTCTTTAGTTTAATTAGATCCCATTTGTCAATTTTGGCTTTTGTTGCCATTGCTTTCGGTGTTTTAGACATCAAGTCTTTGCCCATGCCTATGTCCTGAATGGTATTGCCTAGCTTTTCTTGTAGGGTTTTTATGGTTTTAGGTCTTATGTTTAAG");
	seq_position sp1(0, 735326);
	seq_position sp2(0, 224212723);
	size_t o1 = ref.flatten(sp1);
	size_t o2 = ref.flatten(sp2);

	dna_const_iterator it1 = ref.get_dna(o1-1);
	dna_const_iterator it2 = ref.get_dna(o2);
	SPLOG("%s", assembled.subseq(0, 30).as_string().c_str());
	SPLOG("%s", dna_sequence(it1, it1 + 30).as_string().c_str());
	SPLOG("%s", assembled.subseq(assembled.size() - 30, 30).as_string().c_str());
	SPLOG("%s", dna_sequence(it2 - 30, it2).as_string().c_str());
	std::vector<sv_out> out = call_structural(assembled, dna_slice(it1, it1 + 500), dna_slice(it2 - 500, it2), 500);
	SPLOG("Out size = %d", (int) out.size());
	for(size_t i = 0; i < out.size(); i++) {
		SPLOG("%d: %d) %s (%d", (int) i, 
			(int) (out[i].left_ref - ref.get_dna(0)), 
			dna_sequence(assembled.begin() + out[i].seq_begin, assembled.begin() + out[i].seq_end).as_string().c_str(),
			int (out[i].right_ref - ref.get_dna(0)));
	}
	ASSERT_EQ(out.size(), size_t(1));
	ASSERT_TRUE(out[0].is_structural == false);
}

/*
TEST(call_structural, show)
{
        //                 01234567890123456789012345678901234567890123
        dna_sequence ref1("GACTCACGCTGACTAGCTACCCCATACCAGGGGCATATCATCCA");
        //                 ====*===========***
        dna_sequence seq ("GACTTACGCTGACTAGGGGGCCCATTACGCGGCAGTCCAACAC");
        //                                 ***===============**.\\\\\\\  //
        dna_sequence ref2("CCACAAGCCACAGGATCGAGCCCATTACGCGGCATACCCAACAC");
	std::vector<sv_out> out = call_structural(seq, dna_slice(ref1), dna_slice(ref2), compound_cost(100));
	color_text_buffer ctb;
	pileup p(seq, 150);
	display_var(ctb, seq, out, p);
        file_writer fo(make_path("color_text_buffer.html"));
        fo.print("<html><body>\n");
        ctb.render_as_html(fo);
        fo.print("</body></html>\n");
}
*/


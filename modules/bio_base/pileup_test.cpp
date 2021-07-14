
#include <gtest/gtest.h>
#include "modules/bio_base/pileup.h"

TEST(pileup, test1)
{
	pileup p(dna_sequence("ATTGGATCACTA"), 100);
	ASSERT_TRUE(p.add_read("r1", dna_sequence("GCAT"), "A0AA", true, -1) >= 0);
	ASSERT_TRUE(p.add_read("r2", dna_sequence("GATC"), "AAAA", true, -1) >= 0 );
	ASSERT_TRUE(p.add_read("r3", dna_sequence("AAAA"), "xxxx", false, -1) < 0);
	// 012345
	// ATTGGATCACTA
        //    GCAT
        //     GATC
	ASSERT_EQ(p.depth_at(0), size_t(0));
	ASSERT_EQ(p.depth_at(4), size_t(1));
	ASSERT_EQ(p.depth_at(5), size_t(2));
}


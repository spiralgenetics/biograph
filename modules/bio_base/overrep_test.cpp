
#include "modules/bio_base/overrep.h"
#include <gtest/gtest.h>

TEST(overrep, basic)
{
	overrep_map om(11);
        //                        ACCGTGACCCG
	om.add_overrep(overrep_t("ACCGCTACCCG"_kmer, 100));
	om.add_overrep(overrep_t("ACCGTTACCCG"_kmer, 100));
	om.add_overrep(overrep_t("ACCGTGACCCA"_kmer, 110));
	om.add_overrep(overrep_t("ACCGTGACCCG"_kmer, 120));
	overrep_t o;
	bool r = om.find_near("CGGGTCACGGT"_kmer, o);
	ASSERT_EQ(r, true);
	ASSERT_EQ(o.first, "ACCGTGACCCA"_kmer);
	ASSERT_EQ(o.second, 110);
}


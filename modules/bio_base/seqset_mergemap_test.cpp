#include "modules/bio_base/seqset_mergemap.h"
#include "modules/io/spiral_file.h"
#include "modules/io/spiral_file_mem.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;

TEST(seqset_mergemap, mergemap) {
  spiral_file_mem_storage encoded;
  {
    spiral_file_create_mem c;
    seqset_mergemap_builder b(c.create(), "orig_uuid", "merged_uuid", 20);
    for (size_t i = 0; i < 20; i += 2) {
      b.set(i);
    }
    b.finalize();
    encoded = c.close();
  }

  spiral_file_open_mem opened(encoded);
  seqset_mergemap mergemap(opened.open());

  EXPECT_EQ(mergemap.metadata().orig_seqset_uuid, "orig_uuid");
  EXPECT_EQ(mergemap.metadata().merged_seqset_uuid, "merged_uuid");

  EXPECT_EQ(10, mergemap.get_bitcount().total_bits());
  for (size_t i = 0; i < 20; i++) {
    if (i & 1) {
      EXPECT_FALSE(mergemap.get_bitcount().get(i));
    } else {
      EXPECT_TRUE(mergemap.get_bitcount().get(i));
    }
  }
}

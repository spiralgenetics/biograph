#include "modules/io/membuf.h"
#include "base/base.h"
#include "modules/io/io.h"
#include "modules/io/mmap_buffer.h"

#include <gtest/gtest.h>
#include <sys/mman.h>
#include <random>
#include <unistd.h>

TEST(membuf_test, borrowed_membuf) {
  char x[] = "Borrowed";

  membuf mb(new borrowed_membuf(x, sizeof(x)));
  EXPECT_EQ(mb.size(), sizeof(x));
  EXPECT_EQ(std::string(x), mb.data());
}

TEST(membuf_test, owned_membuf) {
  std::string hi = "Hi!";

  mutable_membuf mb(new owned_membuf(hi.size(), "membuf_test"));
  EXPECT_EQ(hi.size(), mb.size());
  memcpy(mb.mutable_data(), hi.data(), hi.size());
  EXPECT_EQ(0, memcmp(hi.data(), mb.data(), hi.size()));

  membuf not_mutable = mb;
  EXPECT_EQ(hi.size(), not_mutable.size());
  EXPECT_EQ(0, memcmp(hi.data(), not_mutable.data(), hi.size()));
}

TEST(membuf_test, owned_membuf_str) {
  mutable_membuf mb(owned_membuf::from_str("Hello world!", "membuf_test"));
  EXPECT_EQ(mb.str(), "Hello world!");
  EXPECT_EQ(12, mb.size());
}

TEST(membuf_test, subbuf) {
  membuf mb(owned_membuf::from_str("One Two Three", "membuf_test"));
  membuf sub = mb.subbuf(4, 3);
  EXPECT_EQ(3, sub.size());
  EXPECT_EQ("Two", sub.str());
}

TEST(membuf_test, mutablesubbuf) {
  mutable_membuf mb(owned_membuf::from_str("One Two Three", "membuf_test"));
  mutable_membuf sub = mb.subbuf(4, 3);
  EXPECT_EQ(3, sub.size());
  EXPECT_EQ("Two", sub.str());
  memcpy(sub.mutable_data(), "Six", 3);
  EXPECT_EQ("One Six Three", mb.str());
}

// We don't have big files on /share, so it's tough to test caching.
// TODO(nils): Find a way to do this again?
TEST(membuf_test, DISABLED_is_cached) {
  FILE* f = popen("find /share/ -type f -size +20M", "r");
  char buf[1024];
  std::random_device dev;
  std::mt19937 rand_source(dev());
  while (fgets(buf, sizeof(buf), f)) {
    CHECK(*buf);
    buf[strlen(buf) - 1] = 0;

    std::cerr << "Checking if " << buf << " is cached in RAM\n";

    try {
      membuf mb(new mmap_buffer(buf));
      constexpr size_t k_cache_region_size = 20 * 1024 * 1024;
      if (mb.size() > k_cache_region_size) {
        std::uniform_int_distribution<size_t> d(
            0, mb.size() - k_cache_region_size);
        mb = mb.subbuf(d(rand_source), k_cache_region_size);
      }
      if (membuf_cachelist(mb).is_cached_in_memory()) {
        std::cerr << "It is already; skipping\n";
        continue;
      }
      std::cerr << "Attempting to cache " << buf << " in RAM\n";
      membuf_cachelist(mb).cache_in_memory();
      EXPECT_TRUE(membuf_cachelist(mb).is_cached_in_memory());
      pclose(f);
      return;
    } catch (const io_exception& e) {
      std::cerr << e.message() << "\n";
    }
  }
  pclose(f);
  EXPECT_TRUE(false)
      << "Unable to find file that isn't cached already to test\n"
      << "membuf caching.  Try clearing system cache:\n"
      << "  echo 3 > /proc/sys/vm/drop_caches";
}

TEST(membuf_test, discard) {
  constexpr size_t k_mmap_threshold = owned_membuf::k_mmap_threshold;
  constexpr size_t k_tot_size = k_mmap_threshold * 5;
  size_t pagesize = sysconf(_SC_PAGESIZE);

  mutable_membuf mb(new owned_membuf(k_tot_size, "membuf_test"));
  EXPECT_EQ(mb.size(), k_tot_size);
  memset(mb.mutable_data(), 1, k_tot_size);

  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold]);
  EXPECT_EQ(1, mb.data()[2 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[2 * k_mmap_threshold]);
  EXPECT_EQ(1, mb.data()[3 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[3 * k_mmap_threshold]);
  EXPECT_EQ(1, mb.data()[4 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[4 * k_mmap_threshold]);

  // Discard tiny amounts.  We should expect these to do nothing.
  mb.discard_region(mb.mutable_data() + k_mmap_threshold - 1, 1);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold]);

  mb.discard_region(mb.mutable_data() + k_mmap_threshold, 1);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold]);

  mb.discard_region(mb.mutable_data() + k_mmap_threshold - 1, pagesize);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold]);

  mb.discard_region(mb.mutable_data() + k_mmap_threshold + 1,
                       pagesize * 2 - 2);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold + pagesize]);

  // Until it overlaps a whole page.
  mb.discard_region(mb.mutable_data() + k_mmap_threshold + 1,
                       pagesize * 2 - 1);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold + pagesize - 1]);
  EXPECT_EQ(0, mb.data()[1 * k_mmap_threshold + pagesize]);
  EXPECT_EQ(0, mb.data()[1 * k_mmap_threshold + pagesize * 2 - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold + pagesize * 2]);

  // Discard middle section.  We should expect everything in the
  // middle to be gone.
  mb.discard_region(mb.mutable_data() + k_mmap_threshold * 2,
                       k_mmap_threshold);

  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold]);
  EXPECT_EQ(1, mb.data()[2 * k_mmap_threshold - 1]);
  EXPECT_EQ(0, mb.data()[2 * k_mmap_threshold]);
  EXPECT_EQ(0, mb.data()[3 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[3 * k_mmap_threshold]);
  EXPECT_EQ(1, mb.data()[4 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[4 * k_mmap_threshold]);

  // Check for off by 1 errors when discarding near page boundaries.
  mb.discard_region(mb.mutable_data() + k_mmap_threshold + 1,
                       k_mmap_threshold * 3 - 2);

  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold]);
  EXPECT_EQ(0, mb.data()[2 * k_mmap_threshold - 1]);
  EXPECT_EQ(0, mb.data()[2 * k_mmap_threshold]);
  EXPECT_EQ(0, mb.data()[3 * k_mmap_threshold - 1]);
  EXPECT_EQ(0, mb.data()[3 * k_mmap_threshold]);
  EXPECT_EQ(1, mb.data()[4 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[4 * k_mmap_threshold]);

  mb.discard_region(mb.mutable_data() + k_mmap_threshold - 1,
                       k_mmap_threshold * 3 + 2);

  EXPECT_EQ(1, mb.data()[1 * k_mmap_threshold - 1]);
  EXPECT_EQ(0, mb.data()[1 * k_mmap_threshold]);
  EXPECT_EQ(0, mb.data()[2 * k_mmap_threshold - 1]);
  EXPECT_EQ(0, mb.data()[2 * k_mmap_threshold]);
  EXPECT_EQ(0, mb.data()[3 * k_mmap_threshold - 1]);
  EXPECT_EQ(0, mb.data()[3 * k_mmap_threshold]);
  EXPECT_EQ(0, mb.data()[4 * k_mmap_threshold - 1]);
  EXPECT_EQ(1, mb.data()[4 * k_mmap_threshold]);
}

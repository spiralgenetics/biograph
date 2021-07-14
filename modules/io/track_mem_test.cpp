#include "modules/io/track_mem.h"
#include "modules/io/mem_io.h"
#include "modules/io/membuf.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;

class track_mem_test : public Test {
 public:
  void SetUp() override { track_mem::reset_stats(); }

  void TearDown() override {
    for (const auto& type : track_mem::k_types) {
      EXPECT_EQ(type->current_usage(), 0);
      EXPECT_THAT(type->get_detail_usage(), IsEmpty());
    }

    track_mem::reset_stats();
    for (const auto& type : track_mem::k_types) {
      EXPECT_EQ(type->max_usage(), 0);
      EXPECT_EQ(type->current_usage(), 0);
      EXPECT_THAT(type->get_detail_usage(), IsEmpty());
    }

    set_maximum_mem_bytes(48ULL * 1024 * 1024 * 1024);
  }

  static constexpr size_t k_max_untracked_bytes = track_mem::k_max_untracked_bytes;
  static constexpr size_t k_mmap_threshold = owned_membuf::k_mmap_threshold;
};

constexpr size_t track_mem_test::k_max_untracked_bytes;
constexpr size_t track_mem_test::k_mmap_threshold;

TEST_F(track_mem_test, vector_allocator_small_untracked) {
  tracked_vector<char> v(track_alloc("small_vector"));
  v.resize(k_max_untracked_bytes);

  EXPECT_EQ(track_mem::g_malloc_tracker.current_usage(), 0);
  EXPECT_THAT(track_mem::g_malloc_tracker.get_detail_usage(),
              UnorderedElementsAre(Pair("small_vector", 0)));
}

TEST_F(track_mem_test, vector_allocator) {
  tracked_vector<char> v(track_alloc("big_vector"));
  v.resize(k_max_untracked_bytes + 1);

  EXPECT_EQ(track_mem::g_malloc_tracker.current_usage(), k_max_untracked_bytes + 1);
  EXPECT_THAT(track_mem::g_malloc_tracker.get_detail_usage(),
              UnorderedElementsAre(Pair("big_vector", k_max_untracked_bytes + 1)));
}

TEST_F(track_mem_test, small_membuf) {
  membuf m(new owned_membuf(k_mmap_threshold - 1, "small_membuf"));
  EXPECT_EQ(track_mem::g_malloc_tracker.current_usage(), k_mmap_threshold - 1);
  EXPECT_THAT(track_mem::g_malloc_tracker.get_detail_usage(),
              UnorderedElementsAre(Pair("small_membuf", k_mmap_threshold - 1)));
}

TEST_F(track_mem_test, big_membuf) {
  membuf m(new owned_membuf(k_mmap_threshold, "big_membuf"));
  EXPECT_EQ(track_mem::g_malloc_tracker.current_usage(), k_mmap_threshold);
  EXPECT_THAT(track_mem::g_malloc_tracker.get_detail_usage(),
              UnorderedElementsAre(Pair("big_membuf(mmap)", k_mmap_threshold)));
}

TEST_F(track_mem_test, small_mem_io) {
  mem_io m("", track_alloc("small_mem_io"));

  m.reserve(k_max_untracked_bytes / 2);
  EXPECT_EQ(track_mem::g_malloc_tracker.current_usage(), 0);
  EXPECT_THAT(track_mem::g_malloc_tracker.get_detail_usage(),
              UnorderedElementsAre(Pair("small_mem_io", 0)));
}

TEST_F(track_mem_test, big_mem_io) {
  mem_io m("", track_alloc("big_mem_io"));

  m.reserve(k_max_untracked_bytes);
  EXPECT_GE(track_mem::g_malloc_tracker.current_usage(), k_max_untracked_bytes);
  EXPECT_LE(track_mem::g_malloc_tracker.current_usage(), k_max_untracked_bytes * 2);
  EXPECT_THAT(track_mem::g_malloc_tracker.get_detail_usage(),
              UnorderedElementsAre(Pair(
                  "big_mem_io", AllOf(Ge(k_max_untracked_bytes), Le(k_max_untracked_bytes * 2)))));
}

TEST_F(track_mem_test, big_unordered_map) {
  // TODO(nils): Can we rework this test to not require buckets to be the
  // same size as unordered_map::pointer?
  tracked_unordered_map<uint16_t, uint16_t> m(track_alloc("big_unordered_map"));
  size_t elem_size = sizeof(decltype(m)::pointer);
  size_t nelem = (k_max_untracked_bytes + 1) / elem_size + 1;
  m.reserve(nelem * elem_size);
  nelem = m.bucket_count();
  EXPECT_EQ(track_mem::g_malloc_tracker.current_usage(), nelem * elem_size)
      << "Nelem: " << nelem << " Elem size: " << elem_size;
  EXPECT_THAT(track_mem::g_malloc_tracker.get_detail_usage(),
              UnorderedElementsAre(Pair("big_unordered_map", nelem * elem_size)));
}

TEST_F(track_mem_test, does_not_die_under_limit) {
  // This has to be larger than k_min_interesting_hiwat_bytes.
  constexpr size_t k_max_mem = 2ULL * 1024 * 1024 * 1000;

  EXPECT_GT(get_maximum_mem_bytes(), 0);
  set_maximum_mem_bytes(k_max_mem);
  EXPECT_EQ(get_maximum_mem_bytes(), k_max_mem);

  {
    membuf m(new owned_membuf(k_max_mem, "track_mem_limit_test"));
    EXPECT_EQ(m.size(), k_max_mem);
  }
}

class track_mem_death_test : public track_mem_test {};
TEST_F(track_mem_death_test, dies_over_limit) {
  // This has to be larger than k_min_interesting_hiwat_bytes.
  constexpr size_t k_max_mem = 2ULL * 1024 * 1024 * 1000;

  set_maximum_mem_bytes(k_max_mem);

#if !NDEBUG
  // In release mode, it should only warn to the log, not die.
  EXPECT_DEATH(
      {
#endif
        // 1% over triggers death
        membuf m(new owned_membuf(k_max_mem * 102 / 100, "track_mem_death_test"));
#if !NDEBUG
      },
      "exceeds configured maximum");
#endif
}

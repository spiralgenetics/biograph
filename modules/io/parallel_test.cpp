#include "modules/io/parallel.h"
#include "modules/io/make_unique.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace testing;

namespace {

// Returns start + (start+1) + ... + (limit-2) + (limit-1).
size_t sum_from(size_t start, size_t limit) {
  size_t tot = 0;
  // This is faster to do with a formula, but we're more interested in
  // correctness and ease of reading for testing.
  for (size_t i = start; i < limit; ++i) {
    tot += i;
  }
  return tot;
}

}  // namespace

class parallel_test : public testing::Test {
 public:
  void SetUp() override {
    m_orig_threads = get_thread_count();
    m_orig_splits = g_parallel_splits;
  }

  void TearDown() override {
    set_thread_count(m_orig_threads);
    g_parallel_splits = m_orig_splits;
  }

 private:
  int m_orig_threads;
  int m_orig_splits;
};

TEST_F(parallel_test, sum_from) {
  EXPECT_EQ(5, sum_from(5, 6));
  EXPECT_EQ(11, sum_from(5, 7));
  EXPECT_EQ(18, sum_from(5, 8));
}

TEST_F(parallel_test, parallel_for_small) {
  std::mutex mu;
  size_t sum = 0;
  std::map<size_t, size_t> range_size_counts;

  parallel_for(0, 10, [&mu, &sum, &range_size_counts](size_t start, size_t limit) {
    std::lock_guard<std::mutex> l(mu);
    sum += sum_from(start, limit);
    if (!range_size_counts.count(limit - start)) {
      range_size_counts[limit - start] = 0;
    }
    range_size_counts[(limit - start)]++;
  });

  EXPECT_EQ(sum_from(0, 10), sum);
  EXPECT_THAT(range_size_counts, ElementsAre(Pair(1, 10)));
}

TEST_F(parallel_test, parallel_for_big) {
  std::mutex mu;
  size_t sum = 0;
  std::map<size_t, size_t> range_size_counts;
  size_t overall_start = g_parallel_splits / 2;
  size_t overall_limit = g_parallel_splits * 3;

  parallel_for(overall_start, overall_limit,
               [&mu, &sum, &range_size_counts](size_t start, size_t limit) {
                 std::lock_guard<std::mutex> l(mu);
                 sum += sum_from(start, limit);
                 if (!range_size_counts.count(limit - start)) {
                   range_size_counts[limit - start] = 0;
                 }
                 range_size_counts[(limit - start)]++;
               });

  EXPECT_EQ(sum_from(overall_start, overall_limit), sum);
  // It should evenly distribute between 2 elements per chunk and 3
  // elements per chunk.
  EXPECT_THAT(range_size_counts,
              ElementsAre(Pair(2, g_parallel_splits / 2), Pair(3, g_parallel_splits / 2)));
}

namespace {

struct parallel_test_exception : public std::exception {
  const char* what() const noexcept override { return "Test exception"; }
};

}  // namespace

TEST_F(parallel_test, exception_propagation) {
  std::atomic<uint64_t> counter{0};
  unsigned count_to = std::max(1U, std::thread::hardware_concurrency()) * 2 + 1;

  ASSERT_THROW(parallel_for(  //
                   0, std::numeric_limits<uint32_t>::max(),
                   [&](size_t start) {
                     if (counter.fetch_add(1) == count_to) {
                       throw parallel_test_exception();
                     }
                     // Use lots of time; this test should time out if
                     // exceptions aren't propagated immediately.
                     sleep(1);
                   }),
               parallel_test_exception);
  EXPECT_GT(counter.load(), count_to);
}

TEST_F(parallel_test, range_exception_propagation) {
  std::atomic<uint64_t> counter{0};

  ASSERT_GT(g_parallel_splits, 1000)
      << "Sleep time should be increased for this test to be expected.";

  ASSERT_THROW(parallel_for(  //
                   0, std::numeric_limits<uint32_t>::max(),
                   [&](size_t start, size_t limit) {
                     if (counter.fetch_add(1) == 1) {
                       throw parallel_test_exception();
                     }
                     // Use lots of time; this test should time out if
                     // exceptions aren't propagated immediately.
                     sleep(5);
                   }),
               parallel_test_exception);
  EXPECT_GT(counter.load(), 1);
}

TEST_F(parallel_test, progress) {
  std::vector<double> progs;
  constexpr size_t k_num_items = 30;

  set_thread_count(2);

  parallel_for(
      0, k_num_items,
      [](size_t n) { std::this_thread::sleep_for(std::chrono::milliseconds(n * 20)); },
      [&progs](double prog) { progs.push_back(prog); });

  ASSERT_THAT(progs, Not(IsEmpty()));
  EXPECT_THAT(progs[progs.size() - 1], DoubleEq(1)) << PrintToString(progs);
  EXPECT_TRUE(std::is_sorted(progs.begin(), progs.end())) << PrintToString(progs);
  int middle_pos = progs.size() / 2;
  CHECK_LT(middle_pos, progs.size());
  double middle = progs[middle_pos];
  EXPECT_GT(middle, 0.4) << PrintToString(progs) << " at " << middle_pos;
  EXPECT_LT(middle, 0.6) << PrintToString(progs) << " at " << middle_pos;
  EXPECT_GE(progs.size(), k_num_items * 2 / 3) << PrintToString(progs);
  EXPECT_LE(progs.size(), k_num_items * 3 / 2) << PrintToString(progs);
}

TEST_F(parallel_test, subprogress) {
  std::map<int /* difference between consecutive progresses */, unsigned /*count */> pctdiffs;
  constexpr size_t k_num_items = 10;
  constexpr size_t k_num_inside_items = 3;

  int last_pct = 0;

  set_thread_count(2);
  std::string progs = "progs:";

  std::mutex mu;
  size_t prog_updates_expected = 0;
  size_t prog_updates_seen = 0;
  std::condition_variable got_prog_update;

  // Expect a progress update to be delivered shortly
  auto note_prog_expected = [&]() {
    {
      std::unique_lock<std::mutex> l(mu);
      // Wait for previous updates to be receivd by the progress handler.
      got_prog_update.wait(l, [&]() { return prog_updates_expected == prog_updates_seen; });
      ++prog_updates_expected;
    }
  };

  parallel_for(
      0, k_num_items,
      [note_prog_expected](size_t n) {
        // Each of these work items is 10%.

        // Take 0.3 of that 10% (so 3% of the total) and add
        // additional work items.  This should result in 3 jumps of
        // 1% for each of these work items, or 30 jumps of 1% total.
        parallel_for_subprogress(
            0, k_num_inside_items, [note_prog_expected](size_t n2) { note_prog_expected(); },
            0.3);

        // We have .7 left of this work item, or 7% of the total, that
        // we should now add by finishing up this work item.  There
        // are 10 of these work items, so we should get 10 jumps of 7%
        // total.
        note_prog_expected();
      },
      [&](double prog) {
        {
          std::lock_guard<std::mutex> l(mu);
          CHECK_EQ(prog_updates_seen + 1, prog_updates_expected);
          prog_updates_seen = prog_updates_expected;
          got_prog_update.notify_all();
        }

        int pct = prog * 100;
        int diff = pct - last_pct;
        last_pct = pct;
        if (diff) {
          pctdiffs[diff]++;
        }
        progs += " " + std::to_string(prog);
      });

  EXPECT_EQ(100, last_pct) << progs;
  EXPECT_THAT(pctdiffs,
              ElementsAre(Pair(1, k_num_items * k_num_inside_items), Pair(7, k_num_items)))
      << progs;
}

TEST_F(parallel_test, slow_individual) {
  constexpr size_t k_num_items = 64;
  constexpr size_t k_threads = 16;
  std::atomic<size_t> total_sleep{0};

  set_thread_count(k_threads);
  g_parallel_splits = 2;
  auto start = std::chrono::high_resolution_clock::now();
  parallel_for(0, k_num_items, [&total_sleep](size_t n) {
    size_t sleep_ms = n * n;
    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    total_sleep.fetch_add(sleep_ms);
  });
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> elapsed = end - start;

  size_t max_sleep_ms = k_num_items * k_num_items;
  size_t expected_duration = max_sleep_ms * 2.5 + (total_sleep.load() / k_threads);
  EXPECT_LT(elapsed.count(), expected_duration) << "Total expected sleep: " << total_sleep.load()
                                                << " expected duration: " << expected_duration;
}

TEST_F(parallel_test, max_memory) {
  parallel_pool().set_memory_limit(1000);

  std::vector<thread_pool::work_t> worklist;
  std::mutex mu;
  constexpr unsigned k_num_work = 3;
  constexpr unsigned k_num_subwork = 10;

  std::map<int /* stage */, int /* count */> subwork_done;
  size_t mem_used = 0;

  std::vector<std::string> exec_order;

  set_thread_count(k_num_subwork / 2);
  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < k_num_work; ++i) {
    thread_pool::work_t new_work{[&mu, i, start, &mem_used, &exec_order](parallel_state& st) {
      {
        std::lock_guard<std::mutex> l(mu);
        mem_used += 400;
        EXPECT_LE(mem_used, 1000);
        auto so_far = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = so_far - start;
        std::cout << "Starting work " << i << ", elapsed = " << elapsed.count() << "\n";
        exec_order.push_back("start");
      }
      EXPECT_EQ(st.memory_reserved(), 400);
      std::this_thread::sleep_for(std::chrono::milliseconds(2000));
      parallel_for(0, k_num_subwork, [i](size_t j) {
        // Stagger a little bit so we don't get a race condition between finishing work #1 and
        // starting work #2
        std::this_thread::sleep_for(std::chrono::milliseconds(300 + (i * j)));
      });
      {
        std::lock_guard<std::mutex> l(mu);
        CHECK_GE(mem_used, 400);
        mem_used -= 400;
        auto so_far = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = so_far - start;
        std::cout << "Finishing work " << i << ", elapsed = " << elapsed.count() << "\n";
        exec_order.push_back("finish");
      }
    }};
    new_work.reserve_memory = 400;
    worklist.push_back(new_work);
  }
  parallel_pool().execute_worklist(worklist);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsed = end - start;
  size_t expected_time =
      // Critical path should be:
      // sleep(2000ms) for works 0 and 1, runs at once
      2000 +
      // sleep(300ms) * 2 * 2 for subworks 0 and 1 on half as many threads as works
      300 * 2 * 2 +
      // sleep(2000ms) for work 2
      2000 +
      // sleep(300ms) * 2 for subwork 1 on work 2
      300 * 2;

  EXPECT_THAT(elapsed.count(), DoubleNear(expected_time, expected_time * 0.1));
  EXPECT_THAT(exec_order, ElementsAre("start", "start", "finish", "start", "finish", "finish"));
}

struct throws_on_destruct : public parallel_local {
  ~throws_on_destruct() { do_the_throw(); }

  void do_the_throw() {
    if (do_throw_exception) {
      throw(std::runtime_error("Throwing On Destruct"));
    }
  }

  bool do_throw_exception = false;
};

TEST_F(parallel_test, thread_local_exception_DeathTest) {
  set_thread_count(2);
  EXPECT_DEATH(
      {
        parallel_for(0, 5, [](size_t idx, parallel_state& st) {
          throws_on_destruct* tl_var = st.get_local<throws_on_destruct>();
          EXPECT_TRUE(tl_var);
          if (idx == 3) {
            tl_var->do_throw_exception = true;
          }
        });
      },
      "Throwing On Destruct");
}

struct throws_on_flush : public parallel_local {
  void flush() override {
    if (do_throw_exception) {
      throw(std::runtime_error("Throwing On Flush"));
    }
  }

  bool do_throw_exception = false;
};

TEST_F(parallel_test, thread_local_flush_exception) {
  set_thread_count(2);
  ASSERT_THROW(
      {
        parallel_for(0, 5, [](size_t idx, parallel_state& st) {
          throws_on_flush* tl_var = st.get_local<throws_on_flush>();
          EXPECT_TRUE(tl_var);
          if (idx == 3) {
            tl_var->do_throw_exception = true;
          }
        });
      },
      std::runtime_error);
}

#include "modules/io/packed_vector.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/io/log.h"

#include "gtest/gtest.h"

#include <future>
#include <thread>

template <size_t S>
void test_capacity(size_t count, size_t capacity) {
  spiral_file_create_mem c;
  mutable_packed_vector<int, S> pvec(c.create(), count);
  EXPECT_EQ(count, pvec.size());
  EXPECT_EQ(capacity, pvec.capacity());

  spiral_file_open_mem d(c.close());
  packed_vector<int, S> decoded_pvec(d.open());
  EXPECT_EQ(count, decoded_pvec.size());
  EXPECT_EQ(capacity, decoded_pvec.capacity());
}

template <size_t S>
void test_increment() {
  mutable_packed_vector<int, S> pvec(100, "packed_vector_test");
  for (size_t i = 0; i < 100; i++) {
    auto cell = pvec[i];
    // SPLOG("[%zu]: %d", i, (int)cell);
    // pvec.dump();
    EXPECT_EQ(0, (size_t)cell);

    for (size_t j = 0; j < pvec.max_value(); j++) {
      EXPECT_EQ(j, (size_t)cell);
      cell.safe_increment();
      // SPLOG("[%zu][%zu]: %d", i, j, (int)cell);
      // pvec.dump();
    }

    EXPECT_EQ(pvec.max_value(), (size_t)cell);

    for (size_t j = 0; j < 2; j++) {
      cell.safe_increment();
      // SPLOG("[%zu][%zu]: %d", i, j, (int)cell);
      // pvec.dump();
      EXPECT_EQ(pvec.max_value(), (size_t)cell);
    }
  }
}

template <size_t S>
void test_rw() {
  mutable_packed_vector<int, S> pvec(100, "packed_vector_test");
  for (size_t i = 0; i < 100; i++) {
    pvec[i] = i % (pvec.max_value() + 1);
    // SPLOG("[%zu]", i);
    // pvec.dump();
    // SPLOG(" ");
  }

  for (size_t i = 0; i < 100; i++) {
    // SPLOG("[%zu]: %d", i, (int)pvec[i]);
    EXPECT_EQ(i % (pvec.max_value() + 1), pvec[i]);
  }
}

template <typename T, size_t S>
void test_claim() {
  mutable_packed_vector<T, S> pvec(100, "packed_vector_test");
  constexpr size_t k_maxval = ((1ULL) << S) - 1;
  constexpr size_t k_before_maxval = k_maxval - 1;

  for (size_t i = 0; i < 100; ++i) {
    pvec[i] = k_before_maxval;
  }
  for (size_t i = 0; i < 100; ++i) {
    EXPECT_FALSE(pvec[i].safe_increment());
  }

  for (size_t i = 0; i < 100; ++i) {
    EXPECT_EQ(pvec.claim_next_available(i), pvec.size());
  }

  for (size_t i = 0; i < 100; ++i) {
    pvec[i] = k_before_maxval;
    for (size_t j = 0; j < 100; ++j) {
      size_t claimed = pvec.claim_next_available(j);
      if (j <= i) {
        EXPECT_EQ(claimed, i);
        EXPECT_EQ(k_maxval, pvec[i]);
        pvec[i] = k_before_maxval;
      } else {
        EXPECT_EQ(claimed, pvec.size());
        EXPECT_EQ(k_before_maxval, pvec[i]);
      }
    }
    pvec[i] = k_maxval;
  }
}

TEST(packed_vector, basic) {
  test_capacity<1>(1, sizeof(uintmax_t) * (8 / 1));
  test_capacity<1>(2, sizeof(uintmax_t) * (8 / 1));
  test_capacity<1>(4, sizeof(uintmax_t) * (8 / 1));
  test_capacity<1>(8, sizeof(uintmax_t) * (8 / 1));
  test_capacity<1>(32, sizeof(uintmax_t) * (8 / 1));
  test_capacity<1>(64, sizeof(uintmax_t) * (8 / 1));
  test_capacity<1>(65, sizeof(uintmax_t) * (8 / 1) * 2);

  test_capacity<2>(1, sizeof(uintmax_t) * (8 / 2));
  test_capacity<2>(2, sizeof(uintmax_t) * (8 / 2));
  test_capacity<2>(32, sizeof(uintmax_t) * (8 / 2));
  test_capacity<2>(33, sizeof(uintmax_t) * (8 / 2) * 2);

  mutable_packed_vector<int, 2> pvec1(4, "packed_vector_test");
  pvec1[1] = 1;
  pvec1[2] = 2;
  EXPECT_EQ(0, pvec1[0]);
  EXPECT_EQ(1, pvec1[1]);
  EXPECT_EQ(2, pvec1[2]);
  EXPECT_EQ(0, pvec1[3]);
  pvec1[2] = 0;
  EXPECT_EQ(0, pvec1[0]);
  EXPECT_EQ(1, pvec1[1]);
  EXPECT_EQ(0, pvec1[2]);
  EXPECT_EQ(0, pvec1[3]);

  test_rw<2>();
  // test_rw<3>();

  test_increment<2>();
  // test_increment<3>();
  test_increment<4>();

  test_claim<unsigned short, 1>();
  test_claim<unsigned, 1>();
  test_claim<int, 1>();
  test_claim<size_t, 1>();

  test_claim<size_t, 2>();
  test_claim<short, 2>();

  test_claim<size_t, 4>();
  test_claim<unsigned, 4>();

  test_claim<size_t, 8>();
  test_claim<unsigned, 8>();

  test_claim<size_t, 16>();
  test_claim<unsigned, 16>();

  test_claim<size_t, 32>();
  test_claim<unsigned, 32>();
}

TEST(packed_vector, threads) {
  const size_t SIZE = 8000000;
  const size_t NUM_THREADS = 8;

  spiral_file_create_mem c;
  mutable_packed_vector<int, 2> pvec(c.create(), SIZE);
  // pvec.dump();

  std::vector<std::future<void>> futures;
  for (size_t i = 0; i < NUM_THREADS; i++) {
    auto run_thread = [&](size_t i) {
      for (size_t j = 0; j < pvec.size() / NUM_THREADS; j++) {
        auto index = j * NUM_THREADS + i;
        // SPLOG("[%zu][%zu]: %zu", i, j, index);
        pvec[index].safe_increment();
      }
    };
    futures.emplace_back(std::async(std::launch::async, run_thread, i));
  }

  for (auto& future : futures) {
    future.get();
  }

  spiral_file_open_mem d(c.close());
  packed_vector<int, 2> decoded_pvec(d.open());

  // pvec.dump();
  for (size_t j = 0; j < decoded_pvec.size(); j++) {
    EXPECT_EQ(1, (int)decoded_pvec.at(j));
  }
}

TEST(packed_vector, compare_and_swap) {
  const size_t SIZE = 8000000;
  const size_t NUM_THREADS = 8;

  spiral_file_create_mem c;
  mutable_packed_vector<int, 2> pvec(c.create(), SIZE);
  // pvec.dump();

  std::vector<std::future<void>> futures;
  for (size_t i = 0; i < NUM_THREADS; i++) {
    auto run_thread = [&](size_t i) {
      for (size_t j = 0; j < pvec.size() / NUM_THREADS; j++) {
        auto index = j * NUM_THREADS + i;
        // SPLOG("[%zu][%zu]: %zu", i, j, index);
        for (;;) {
          int old_value = pvec.at(index);
          int new_value = old_value + 1;
          if (pvec.at(index).compare_and_swap(old_value, new_value)) {
            break;
          }
        }
      }
    };
    futures.emplace_back(std::async(std::launch::async, run_thread, i));
  }

  for (auto& future : futures) {
    future.get();
  }

  spiral_file_open_mem d(c.close());
  packed_vector<int, 2> decoded_pvec(d.open());

  // pvec.dump();
  for (size_t j = 0; j < decoded_pvec.size(); j++) {
    EXPECT_EQ(1, (int)decoded_pvec.at(j));
  }
}


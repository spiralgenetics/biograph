#include "benchmark/benchmark.h"
#include "modules/io/membuf.h"
#include "modules/io/parallel.h"
#include "modules/io/progress.h"
#include "modules/io/utils.h"

#include <iomanip>
#include <iostream>
#include <random>

#ifdef _OPENMP

#include <parallel/algorithm>

// Benchmark a bunch of GNU parallel sort algorithms.

namespace {

constexpr size_t k_MB = 1024 * 1024;

// Maximum amount of RAM we should use.
constexpr size_t k_max_ram_gb = 230;

// Use up to full RAM for qsort, since it sorts in place.
constexpr size_t k_max_qsort_mb = k_max_ram_gb * 1024;
constexpr size_t k_min_qsort_mb = k_max_qsort_mb / 10;

// Non-nested qsorts take a long time, so don't run them on big dataets.
constexpr size_t k_max_slow_qsort_mb = 20 * 1024;
constexpr size_t k_min_slow_qsort_mb = k_max_slow_qsort_mb / 4;

// Merge needs to copy all the data, so we can only use half the RAM.
constexpr size_t k_max_merge_mb = k_max_ram_gb * 1024 / 2;
constexpr size_t k_min_merge_mb = k_max_merge_mb / 10;

// Do range steps in factors of 3; this makes it less likely we'll get
// any benefit by hitting sizes that are factors of 2.
constexpr unsigned k_range_multiplier = 3;

mutable_membuf elem_buf;
uint64_t* elems_begin = nullptr;
uint64_t* elems_end = nullptr;
size_t elems_size = 0;

void init_elems(size_t size) {
  elem_buf = mutable_membuf();
  size_t new_bytes = size * k_MB;
  CHECK_LT(new_bytes, get_system_mem());
  elem_buf = mutable_membuf(new owned_membuf(size * k_MB, "parallel_benchmark"));
  elems_size = elem_buf.size() / sizeof(uint64_t);

  elems_begin = reinterpret_cast<uint64_t*>(elem_buf.mutable_data());
  elems_end = elems_begin + elems_size;

  parallel_for(  //
      0, elems_size, [&](size_t start, size_t limit) {
        std::random_device rand_dev;
        std::mt19937_64 random_source(rand_dev());

        auto end = elems_begin + limit;

        for (auto it = elems_begin + start; it != end; ++it) {
          *it = random_source();
        }
      });
}

}  // namespace

static void BM_parallel_mergesort_exact(benchmark::State& state) {
  while (state.KeepRunning()) {
    state.PauseTiming();
    init_elems(state.range(0));
    state.ResumeTiming();

    __gnu_parallel::sort(elems_begin, elems_end, __gnu_parallel::multiway_mergesort_exact_tag());
  }
  state.SetBytesProcessed(elems_size * sizeof(uint64_t) * state.iterations());
}

BENCHMARK(BM_parallel_mergesort_exact)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(k_range_multiplier)
    ->Range(k_min_merge_mb, k_max_merge_mb);

static void BM_parallel_mergesort_sampling(benchmark::State& state) {
  while (state.KeepRunning()) {
    state.PauseTiming();
    init_elems(state.range(0));
    state.ResumeTiming();

    __gnu_parallel::sort(elems_begin, elems_end, __gnu_parallel::multiway_mergesort_sampling_tag());
  }
  state.SetBytesProcessed(elems_size * sizeof(uint64_t) * state.iterations());
}

BENCHMARK(BM_parallel_mergesort_sampling)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(k_range_multiplier)
    ->Range(k_min_merge_mb, k_max_merge_mb);

static void BM_parallel_qsort_nested(benchmark::State& state) {
  while (state.KeepRunning()) {
    state.PauseTiming();
    init_elems(state.range(0));
    state.ResumeTiming();

    omp_set_nested(1);
    __gnu_parallel::sort(elems_begin, elems_end, __gnu_parallel::quicksort_tag());
  }
  state.SetBytesProcessed(elems_size * sizeof(uint64_t) * state.iterations());
}

BENCHMARK(BM_parallel_qsort_nested)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(k_range_multiplier)
    ->Range(k_min_qsort_mb, k_max_qsort_mb);

static void BM_parallel_qsort_balanced_nested(benchmark::State& state) {
  while (state.KeepRunning()) {
    state.PauseTiming();
    init_elems(state.range(0));
    state.ResumeTiming();

    omp_set_nested(1);
    __gnu_parallel::sort(elems_begin, elems_end, __gnu_parallel::balanced_quicksort_tag());
  }
  state.SetBytesProcessed(elems_size * sizeof(uint64_t) * state.iterations());
}

BENCHMARK(BM_parallel_qsort_balanced_nested)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(k_range_multiplier)
    ->Range(k_min_qsort_mb, k_max_qsort_mb);

static void BM_parallel_qsort(benchmark::State& state) {
  while (state.KeepRunning()) {
    state.PauseTiming();
    init_elems(state.range(0));
    state.ResumeTiming();

    omp_set_nested(0);
    __gnu_parallel::sort(elems_begin, elems_end, __gnu_parallel::quicksort_tag());
  }
  state.SetBytesProcessed(elems_size * sizeof(uint64_t) * state.iterations());
}

BENCHMARK(BM_parallel_qsort)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(k_range_multiplier)
    ->Range(k_min_slow_qsort_mb, k_max_slow_qsort_mb);

static void BM_parallel_qsort_balanced(benchmark::State& state) {
  while (state.KeepRunning()) {
    state.PauseTiming();
    init_elems(state.range(0));
    state.ResumeTiming();

    omp_set_nested(0);
    __gnu_parallel::sort(elems_begin, elems_end, __gnu_parallel::balanced_quicksort_tag());
  }
  state.SetBytesProcessed(elems_size * sizeof(uint64_t) * state.iterations());
}

BENCHMARK(BM_parallel_qsort_balanced)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(k_range_multiplier)
    ->Range(k_min_slow_qsort_mb, k_max_slow_qsort_mb);

#endif

BENCHMARK_MAIN();
